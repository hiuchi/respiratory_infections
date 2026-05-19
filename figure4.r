library(DESeq2)
library(patchwork)
library(scales)
library(tidyverse)

figure_palettes <- list(
  infection = c(
    Mock = "#8A8A8A",
    Flu = "#3C78B5",
    COVID = "#D95A4E"
  ),
  age = c(
    Young = "#2FA37A",
    Aged = "#7B62B3"
  )
)

figure_strip_fill <- "#F4F4F4"
figure_strip_border <- "#D0D0D0"

parse_time_hours <- function(time_label) {
  out <- rep(Inf, length(time_label))
  has_value <- !is.na(time_label)
  is_hours <- has_value & str_detect(time_label, "^[0-9]+h$")
  is_days <- has_value & str_detect(time_label, "^[0-9]+d$")
  out[is_hours] <- as.numeric(str_remove(time_label[is_hours], "h$"))
  out[is_days] <- 24 * as.numeric(str_remove(time_label[is_days], "d$"))
  out
}

parse_sample_metadata <- function(samples) {
  sample_ids <- as.character(samples)
  match <- str_match(
    sample_ids,
    "^([^-]+)-([^-]+)-([0-9]+[hd])(?:-(.+))?$"
  )

  tibble(
    sample_id = sample_ids,
    series_code = match[, 2],
    condition = match[, 3],
    time_label = match[, 4],
    replicate = match[, 5]
  ) |>
    mutate(
      series = case_when(
        series_code == "y" ~ "Young",
        series_code == "a" ~ "Aged",
        !is.na(series_code) ~ str_to_title(series_code),
        TRUE ~ NA_character_
      ),
      time_hours = parse_time_hours(time_label),
      replicate_num = parse_number(replicate),
      is_valid = !is.na(series_code) &
        !is.na(condition) &
        !is.na(time_label) &
        is.finite(time_hours)
    )
}

build_figure_output <- function(project_dir, filename) {
  output_dir <- file.path(project_dir, "output", "plots")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  file.path(output_dir, filename)
}

save_figure_pdf <- function(plot_obj, output_pdf, figure_width, figure_height) {
  ggsave(
    filename = output_pdf,
    plot = plot_obj,
    width = figure_width,
    height = figure_height,
    units = "in",
    limitsize = FALSE
  )
}

tpm_file <- "/Users/hiuchi/Dropbox/Research/covid/260305_analysis/res/star_salmon/salmon.merged.gene_tpm.tsv"
input_dir <- "/Users/hiuchi/Dropbox/Research/covid/260305_analysis/res/star_salmon/deseq2_qc"
project_dir <- "/Users/hiuchi/Dropbox/Research/covid/260309_analysis"
target_time <- "2d"
max_padj <- 0.05
min_base_mean <- 20
top_heatmap_n <- list(
  virus = c(same = 30, specific = 30),
  age = c(same = 30, specific = 30)
)
top_panel_n <- list(
  virus = c(same = 3, specific = 3),
  age = c(same = 3, specific = 3)
)
panel_priority_genes <- list(
  virus = list(
    same = c("S100a8", "S100a9", "Ccl20", "Csf3", "Osm", "Il1r2", "Il36g", "Nlrp12"),
    specific = character()
  ),
  age = list(
    same = c("Igha", "Igkc", "Cfh"),
    specific = character()
  )
)

output_pdf <- Sys.getenv("FIGURE4_OUTPUT_PDF", unset = "")
if (output_pdf == "") {
  output_pdf <- build_figure_output(project_dir, "figure4.pdf")
}

figure_width <- 14
figure_height <- 12.5

plot_title_size <- 12.5
axis_title_size <- 11
axis_text_size <- 10
legend_text_size <- 9.5
strip_text_size <- 10.8
heatmap_x_text_size <- 9
heatmap_y_text_size <- 7.6
heatmap_tile_width <- 0.95

infection_levels <- c("Mock", "Flu", "COVID")
age_levels <- c("Young", "Aged")

condition_palette <- figure_palettes$infection
age_palette <- figure_palettes$age
heatmap_palette <- c(
  low = "#4E79A7",
  mid = "#F5F2EC",
  high = "#C94F4F"
)

module_titles <- c(
  virus = "Flu-COVID difference genes",
  age = "Young-Aged difference genes"
)

pattern_labels <- c(
  same = "Shared",
  specific = "Specific"
)

parse_sample_info <- function(samples, infection_levels, age_levels) {
  sample_info <- parse_sample_metadata(samples) |>
    filter(
      is_valid,
      series %in% age_levels,
      condition %in% infection_levels,
      !is.na(replicate_num)
    )

  time_levels <- sample_info |>
    distinct(time_label, time_hours) |>
    arrange(time_hours) |>
    pull(time_label)

  sample_info |>
    transmute(
      sample = sample_id,
      age = factor(series, levels = age_levels),
      infection = factor(condition, levels = infection_levels),
      time = factor(time_label, levels = time_levels),
      time_hours
    )
}

make_gene_lookup <- function(gene_sets) {
  purrr::imap_dfr(gene_sets, function(pattern_sets, module_name) {
    purrr::imap_dfr(pattern_sets, function(genes, pattern_name) {
      tibble(module = module_name, pattern = pattern_name, gene_name = genes)
    })
  })
}

add_module_metadata <- function(df, module_lookup, module_levels, gene_levels, pattern_labels) {
  df |>
    left_join(module_lookup, by = "gene_name", relationship = "many-to-many") |>
    filter(!is.na(module)) |>
    mutate(
      gene_name = factor(gene_name, levels = unique(gene_levels)),
      module = factor(module, levels = module_levels),
      pattern = factor(pattern, levels = names(pattern_labels), labels = pattern_labels)
    )
}

drop_extra_assays <- function(dds_subset) {
  SummarizedExperiment::assays(dds_subset) <- S4Vectors::SimpleList(
    counts = counts(dds_subset)
  )
  dds_subset
}

compute_base_mean <- function(dds, sample_ids, label) {
  normalized_counts <- counts(dds[, sample_ids], normalized = TRUE)

  tibble(
    gene_id = rownames(normalized_counts),
    !!paste0("baseMean_", label) := rowMeans(normalized_counts)
  )
}

apply_base_mean_filter <- function(df, base_mean_cols, min_base_mean) {
  keep <- rep(TRUE, nrow(df))

  for (col_name in base_mean_cols) {
    keep <- keep & !is.na(df[[col_name]]) & df[[col_name]] >= min_base_mean
  }

  df[keep, , drop = FALSE]
}

run_contrast <- function(dds, coldata, design_formula, coef, result_label) {
  dds_subset <- dds[, coldata$sample]
  dds_subset <- drop_extra_assays(dds_subset)

  coldata_df <- coldata |>
    select(-sample) |>
    as.data.frame()
  rownames(coldata_df) <- coldata$sample

  colData(dds_subset) <- S4Vectors::DataFrame(coldata_df)
  design(dds_subset) <- design_formula
  dds_subset <- DESeq(dds_subset, quiet = TRUE)

  as.data.frame(results(dds_subset, name = coef)) |>
    rownames_to_column("gene_id") |>
    as_tibble() |>
    transmute(
      gene_id,
      !!paste0("log2FoldChange_", result_label) := log2FoldChange,
      !!paste0("padj_", result_label) := padj
    )
}

build_infection_coldata <- function(sample_info, target_time, age_value) {
  sample_info |>
    filter(
      time == target_time,
      age == age_value,
      infection %in% c("Flu", "COVID")
    ) |>
    transmute(
      sample,
      infection = factor(as.character(infection), levels = c("Flu", "COVID"))
    )
}

build_age_coldata <- function(sample_info, target_time, infection_value) {
  sample_info |>
    filter(
      time == target_time,
      infection == infection_value,
      age %in% c("Young", "Aged")
    ) |>
    transmute(
      sample,
      age = factor(as.character(age), levels = c("Young", "Aged"))
    )
}

rank_paired_deg_hits <- function(
  result_df,
  gene_lookup,
  module_name,
  label_a,
  label_b,
  base_mean_cols,
  max_padj,
  min_base_mean
) {
  lfc_a <- paste0("log2FoldChange_", label_a)
  lfc_b <- paste0("log2FoldChange_", label_b)
  padj_a <- paste0("padj_", label_a)
  padj_b <- paste0("padj_", label_b)

  result_df |>
    left_join(gene_lookup, by = "gene_id") |>
    mutate(
      gene_name = if_else(is.na(gene_name) | gene_name == "", gene_id, gene_name),
      shared_abs_lfc = pmin(abs(.data[[lfc_a]]), abs(.data[[lfc_b]])),
      mean_abs_lfc = (abs(.data[[lfc_a]]) + abs(.data[[lfc_b]])) / 2,
      shared_padj = pmax(.data[[padj_a]], .data[[padj_b]], na.rm = TRUE),
      specificity_score = abs(abs(.data[[lfc_a]]) - abs(.data[[lfc_b]])),
      both_significant = .data[[padj_a]] < max_padj & .data[[padj_b]] < max_padj,
      either_significant = .data[[padj_a]] < max_padj | .data[[padj_b]] < max_padj,
      same_hit = both_significant & sign(.data[[lfc_a]]) == sign(.data[[lfc_b]]),
      specific_hit = either_significant & !same_hit,
      pattern = case_when(
        same_hit ~ "same",
        specific_hit ~ "specific",
        TRUE ~ NA_character_
      ),
      rank_primary = if_else(
        pattern == "specific",
        specificity_score,
        shared_abs_lfc
      ),
      rank_secondary = if_else(
        pattern == "specific",
        0,
        mean_abs_lfc
      ),
      rank_padj = if_else(
        pattern == "specific",
        0,
        shared_padj
      )
    ) |>
    filter(
      !is.na(pattern),
      !is.na(.data[[padj_a]]),
      !is.na(.data[[padj_b]]),
      sign(.data[[lfc_a]]) != 0,
      sign(.data[[lfc_b]]) != 0
    ) |>
    apply_base_mean_filter(
      base_mean_cols = base_mean_cols,
      min_base_mean = min_base_mean
    ) |>
    arrange(
      pattern,
      desc(rank_primary),
      desc(rank_secondary),
      rank_padj,
      gene_name
    ) |>
    distinct(gene_name, .keep_all = TRUE) |>
    mutate(module = module_name)
}

select_top_genes <- function(candidate_table, pattern_name, n_genes, module_name, use_name) {
  genes <- candidate_table |>
    filter(pattern == pattern_name) |>
    slice_head(n = n_genes) |>
    pull(gene_name)

  if (length(genes) == 0) {
    stop("No genes found for module ", module_name, " and use ", use_name, ".")
  }

  genes
}

select_panel_genes <- function(candidate_table, pattern_name, priority_genes, n_genes, module_name) {
  candidate_genes <- candidate_table |>
    filter(pattern == pattern_name) |>
    pull(gene_name)

  genes <- priority_genes[priority_genes %in% candidate_genes]
  if (length(genes) < n_genes) {
    genes <- c(genes, setdiff(candidate_genes, genes))
  }

  if (length(genes) < n_genes) {
    stop("Not enough panel genes found for module ", module_name, " and pattern ", pattern_name, ".")
  }

  genes[seq_len(n_genes)]
}

select_gene_sets_by_pattern <- function(candidate_table, n_by_pattern, module_name, use_name) {
  purrr::imap(n_by_pattern, function(n_genes, pattern_name) {
    select_top_genes(candidate_table, pattern_name, n_genes, module_name, use_name)
  })
}

select_panel_sets_by_pattern <- function(candidate_table, n_by_pattern, priority_by_pattern, module_name) {
  purrr::imap(n_by_pattern, function(n_genes, pattern_name) {
    select_panel_genes(
      candidate_table = candidate_table,
      pattern_name = pattern_name,
      priority_genes = priority_by_pattern[[pattern_name]],
      n_genes = n_genes,
      module_name = module_name
    )
  })
}

select_figure4_gene_sets <- function(
  dds,
  sample_info,
  gene_lookup,
  target_time,
  max_padj,
  min_base_mean,
  top_heatmap_n,
  top_panel_n,
  panel_priority_genes
) {
  infected_2d_base_means <- list(
    compute_base_mean(
      dds,
      sample_info |>
        filter(infection == "Flu", time == target_time) |>
        pull(sample),
      "Flu"
    ),
    compute_base_mean(
      dds,
      sample_info |>
        filter(infection == "COVID", time == target_time) |>
        pull(sample),
      "COVID"
    )
  ) |>
    purrr::reduce(inner_join, by = "gene_id")

  virus_candidates <- list(
    run_contrast(
      dds = dds,
      coldata = build_infection_coldata(sample_info, target_time, "Young"),
      design_formula = ~ infection,
      coef = "infection_COVID_vs_Flu",
      result_label = "virus_young"
    ),
    run_contrast(
      dds = dds,
      coldata = build_infection_coldata(sample_info, target_time, "Aged"),
      design_formula = ~ infection,
      coef = "infection_COVID_vs_Flu",
      result_label = "virus_aged"
    )
  ) |>
    purrr::reduce(inner_join, by = "gene_id") |>
    left_join(infected_2d_base_means, by = "gene_id") |>
    rank_paired_deg_hits(
      gene_lookup = gene_lookup,
      module_name = "virus",
      label_a = "virus_young",
      label_b = "virus_aged",
      base_mean_cols = c("baseMean_Flu", "baseMean_COVID"),
      max_padj = max_padj,
      min_base_mean = min_base_mean
    )

  age_candidates <- list(
    run_contrast(
      dds = dds,
      coldata = build_age_coldata(sample_info, target_time, "Flu"),
      design_formula = ~ age,
      coef = "age_Aged_vs_Young",
      result_label = "age_flu"
    ),
    run_contrast(
      dds = dds,
      coldata = build_age_coldata(sample_info, target_time, "COVID"),
      design_formula = ~ age,
      coef = "age_Aged_vs_Young",
      result_label = "age_covid"
    )
  ) |>
    purrr::reduce(inner_join, by = "gene_id") |>
    left_join(infected_2d_base_means, by = "gene_id") |>
    rank_paired_deg_hits(
      gene_lookup = gene_lookup,
      module_name = "age",
      label_a = "age_flu",
      label_b = "age_covid",
      base_mean_cols = c("baseMean_Flu", "baseMean_COVID"),
      max_padj = max_padj,
      min_base_mean = min_base_mean
    )

  heatmap_gene_sets <- list(
    virus = select_gene_sets_by_pattern(virus_candidates, top_heatmap_n[["virus"]], "virus", "heatmap"),
    age = select_gene_sets_by_pattern(age_candidates, top_heatmap_n[["age"]], "age", "heatmap")
  )

  list(
    heatmap = heatmap_gene_sets,
    panel = list(
      virus = select_panel_sets_by_pattern(
        candidate_table = virus_candidates,
        n_by_pattern = top_panel_n[["virus"]],
        priority_by_pattern = panel_priority_genes$virus,
        "virus"
      ),
      age = select_panel_sets_by_pattern(
        candidate_table = age_candidates,
        n_by_pattern = top_panel_n[["age"]],
        priority_by_pattern = panel_priority_genes$age,
        "age"
      )
    )
  )
}

summarise_mean_expr <- function(df, group_vars, value_col = "log2_tpm") {
  df |>
    group_by(across(all_of(group_vars))) |>
    summarise(mean_expr = mean(.data[[value_col]]), .groups = "drop")
}

module_plot_theme <- function(
  axis_title_size,
  axis_text_size,
  legend_text_size,
  plot_title_size,
  strip_text_size
) {
  theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(size = legend_text_size),
      strip.placement = "outside",
      strip.switch.pad.grid = grid::unit(0, "pt"),
      strip.switch.pad.wrap = grid::unit(0, "pt"),
      strip.background = element_rect(fill = figure_strip_fill, colour = figure_strip_border),
      strip.text = element_text(face = "bold", size = strip_text_size),
      strip.text.y.left = element_text(
        angle = 90,
        face = "bold",
        size = strip_text_size * 0.82,
        lineheight = 0.9
      ),
      plot.title = element_text(face = "bold", size = plot_title_size),
      plot.margin = margin(3, 6, 3, 6),
      axis.text.x = element_text(size = axis_text_size),
      axis.text.y = element_text(size = axis_text_size)
    )
}

build_ratio_heatmap <- function(
  module_df,
  panel_var,
  panel_levels,
  ratio_var,
  numerator_level,
  denominator_level
) {
  ratio_df <- module_df |>
    transmute(
      pattern,
      gene_name,
      panel_group = factor(as.character(.data[[panel_var]]), levels = panel_levels),
      ratio_group = as.character(.data[[ratio_var]]),
      time,
      time_hours,
      mean_expr
    ) |>
    filter(!is.na(panel_group)) |>
    filter(ratio_group %in% c(numerator_level, denominator_level)) |>
    pivot_wider(names_from = ratio_group, values_from = mean_expr) |>
    filter(!is.na(.data[[numerator_level]]), !is.na(.data[[denominator_level]])) |>
    mutate(
      ratio_value = log2(
        (.data[[numerator_level]] + 1) / (.data[[denominator_level]] + 1)
      )
    )

  columns <- ratio_df |>
    distinct(panel_group, time, time_hours) |>
    mutate(
      column_id = paste(panel_group, time, sep = "__"),
      column_label = as.character(time)
    ) |>
    arrange(panel_group, time_hours)

  ratio_df |>
    left_join(columns, by = c("panel_group", "time", "time_hours")) |>
    mutate(
      column_id = factor(column_id, levels = columns$column_id),
      panel_group = factor(panel_group, levels = levels(columns$panel_group)),
      pattern = factor(pattern, levels = unique(pattern))
    )
}

make_heatmap_plot <- function(
  heatmap_df,
  module_levels,
  module_title = NULL,
  heatmap_palette,
  y_label_margin_right,
  plot_title_size,
  strip_text_size,
  heatmap_x_text_size,
  heatmap_y_text_size,
  legend_text_size
) {
  heatmap_labels <- heatmap_df |>
    distinct(column_id, column_label) |>
    arrange(column_id)
  heatmap_label_lookup <- setNames(heatmap_labels$column_label, heatmap_labels$column_id)

  module_heatmap_df <- heatmap_df |>
    mutate(gene_name = factor(gene_name, levels = rev(module_levels)))

  ggplot(module_heatmap_df, aes(column_id, gene_name, fill = ratio_value)) +
    geom_tile(
      width = heatmap_tile_width,
      height = 0.95,
      colour = "white",
      linewidth = 0.25
    ) +
    facet_grid(
      rows = vars(pattern),
      cols = vars(panel_group),
      scales = "free",
      space = "free"
    ) +
    scale_x_discrete(labels = heatmap_label_lookup, expand = expansion(mult = c(0, 0))) +
    scale_y_discrete(expand = expansion(mult = c(0, 0))) +
    scale_fill_gradient2(
      low = heatmap_palette[["low"]],
      mid = heatmap_palette[["mid"]],
      high = heatmap_palette[["high"]],
      midpoint = 0,
      limits = c(-2.5, 2.5),
      oob = squish,
      name = "log2\nratio"
    ) +
    labs(title = module_title, x = NULL, y = NULL) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      strip.background = element_rect(fill = figure_strip_fill, colour = figure_strip_border),
      strip.text = element_text(face = "bold", size = strip_text_size),
      strip.text.y = element_text(angle = -90),
      axis.text.x = element_text(
        size = heatmap_x_text_size,
        angle = 0,
        vjust = 1,
        hjust = 0.5
      ),
      axis.text.y = element_text(
        size = heatmap_y_text_size,
        hjust = 1,
        vjust = 0.5,
        margin = margin(r = y_label_margin_right)
      ),
      axis.ticks.y = element_blank(),
      legend.title = element_text(size = legend_text_size),
      legend.text = element_text(size = legend_text_size),
      plot.title = element_text(face = "bold", size = plot_title_size),
      plot.margin = margin(4, 6, 2, 9)
    )
}

make_heatmap_scale_legend_plot <- function(
  heatmap_palette,
  legend_text_size
) {
  legend_limits <- c(-2.5, 2.5)
  legend_bar_limits <- c(-0.5, 0.5)
  legend_breaks <- c(-2, -1, 0, 1, 2)
  legend_df <- tibble(
    z = seq(legend_limits[1], legend_limits[2], length.out = 400),
    x = scales::rescale(z, to = legend_bar_limits, from = legend_limits),
    y = 1
  )
  axis_breaks <- scales::rescale(
    legend_breaks,
    to = legend_bar_limits,
    from = legend_limits
  )

  ggplot(legend_df, aes(x, y, fill = z)) +
    geom_tile(height = 1.8) +
    scale_fill_gradient2(
      low = heatmap_palette[["low"]],
      mid = heatmap_palette[["mid"]],
      high = heatmap_palette[["high"]],
      midpoint = 0,
      limits = legend_limits,
      guide = "none"
    ) +
    scale_x_continuous(
      breaks = axis_breaks,
      labels = legend_breaks,
      limits = legend_limits,
      expand = expansion(mult = c(0, 0))
    ) +
    labs(x = "log2 ratio", y = NULL) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(size = legend_text_size),
      axis.title.x = element_text(size = legend_text_size, face = "bold"),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      plot.margin = margin(2, 10, 2, 4)
    )
}

make_module_plot <- function(
  module_name,
  expr_long,
  module_levels,
  module_title,
  condition_palette,
  age_palette,
  plot_time_levels,
  axis_title_size,
  axis_text_size,
  legend_text_size,
  plot_title_size,
  strip_text_size
) {
  base_theme <- module_plot_theme(
    axis_title_size = axis_title_size,
    axis_text_size = axis_text_size,
    legend_text_size = legend_text_size,
    plot_title_size = plot_title_size,
    strip_text_size = strip_text_size
  )

  module_data <- expr_long |>
    filter(module == module_name) |>
    mutate(
      gene_name = factor(gene_name, levels = module_levels),
      age = factor(age, levels = c("Young", "Aged")),
      pattern = factor(pattern, levels = pattern_labels),
      row_label = paste(as.character(pattern), gene_name, sep = "\n")
    )
  row_label_levels <- module_data |>
    distinct(gene_name, row_label) |>
    arrange(gene_name) |>
    pull(row_label)
  module_data <- module_data |>
    mutate(row_label = factor(row_label, levels = row_label_levels))

  if (module_name == "virus") {
    module_means <- summarise_mean_expr(
      module_data,
      c("row_label", "gene_name", "age", "infection", "time", "time_hours")
    )

    return(
      ggplot(module_data, aes(time, log2_tpm, colour = infection)) +
        geom_point(
          position = position_jitter(width = 0.08, height = 0),
          size = 1.15,
          alpha = 0.22
        ) +
        geom_line(
          data = module_means |>
            filter(infection != "Mock"),
          aes(y = mean_expr, group = infection),
          linewidth = 0.65,
          alpha = 0.95
        ) +
        geom_point(
          data = module_means,
          aes(y = mean_expr),
          size = 1.65,
          alpha = 0.95
        ) +
        facet_grid(
          rows = vars(row_label),
          cols = vars(age),
          scales = "free_y",
          switch = "y"
        ) +
        scale_colour_manual(values = condition_palette) +
        scale_x_discrete(limits = plot_time_levels, drop = FALSE) +
        expand_limits(y = 0) +
        scale_y_continuous(
          expand = expansion(mult = c(0, 0.05)),
          sec.axis = dup_axis(name = "log2(TPM + 1)")
        ) +
        guides(
          colour = guide_legend(
            override.aes = list(alpha = 1, size = 2.6)
          )
        ) +
        labs(
          title = module_title,
          x = NULL,
          y = NULL
        ) +
        base_theme +
        theme(
          legend.position = "top",
          axis.text.y.left = element_blank(),
          axis.ticks.y.left = element_blank(),
          axis.line.y.left = element_blank(),
          axis.text.y.right = element_text(size = axis_text_size),
          axis.ticks.y.right = element_line(),
          axis.title.y.right = element_text(
            size = axis_title_size,
            angle = -90,
            margin = margin(l = 6)
          )
        )
    )
  }

  module_data <- module_data |>
    mutate(infection = factor(infection, levels = c("Mock", "Flu", "COVID")))

  module_means <- summarise_mean_expr(
    module_data,
    c("row_label", "gene_name", "infection", "age", "time", "time_hours")
  )
  line_means <- module_means |>
    add_count(gene_name, infection, age, name = "n_points") |>
    filter(n_points > 1)

  ggplot(module_data, aes(time, log2_tpm, colour = age)) +
    geom_point(
      position = position_jitter(width = 0.08, height = 0),
      size = 1.15,
      alpha = 0.22
    ) +
    geom_line(
      data = line_means,
      aes(y = mean_expr, group = age),
      linewidth = 0.65,
      alpha = 0.95
    ) +
    geom_point(
      data = module_means,
      aes(y = mean_expr),
      size = 1.65,
      alpha = 0.95
    ) +
    facet_grid(
      rows = vars(row_label),
      cols = vars(infection),
      scales = "free_y",
      switch = "y"
    ) +
    scale_colour_manual(values = age_palette) +
    scale_x_discrete(limits = plot_time_levels, drop = FALSE) +
    expand_limits(y = 0) +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.05)),
      sec.axis = dup_axis(name = "log2(TPM + 1)")
    ) +
    guides(
      colour = guide_legend(
        override.aes = list(alpha = 1, size = 2.6)
      )
    ) +
    labs(
      title = module_title,
      x = NULL,
      y = NULL
    ) +
    base_theme +
    theme(
      legend.position = "top",
      axis.text.y.left = element_blank(),
      axis.ticks.y.left = element_blank(),
      axis.line.y.left = element_blank(),
      axis.text.y.right = element_text(size = axis_text_size),
      axis.ticks.y.right = element_line(),
      axis.title.y.right = element_text(
        size = axis_title_size,
        angle = -90,
        margin = margin(l = 6)
      )
    )
}

tpm_table <- readr::read_tsv(tpm_file, show_col_types = FALSE)
raw_sample_columns <- setdiff(names(tpm_table), c("gene_id", "gene_name"))
sample_info <- parse_sample_info(raw_sample_columns, infection_levels, age_levels)
sample_columns <- sample_info$sample
gene_lookup <- tpm_table |>
  distinct(gene_id, gene_name) |>
  mutate(
    gene_name = if_else(is.na(gene_name) | gene_name == "", gene_id, gene_name)
  )

load(file.path(input_dir, "deseq2.dds.RData"))

gene_sets <- select_figure4_gene_sets(
  dds = dds,
  sample_info = sample_info,
  gene_lookup = gene_lookup,
  target_time = target_time,
  max_padj = max_padj,
  min_base_mean = min_base_mean,
  top_heatmap_n = top_heatmap_n,
  top_panel_n = top_panel_n,
  panel_priority_genes = panel_priority_genes
)

heatmap_gene_sets <- gene_sets$heatmap
panel_gene_sets <- gene_sets$panel
heatmap_module_lookup <- make_gene_lookup(heatmap_gene_sets)
panel_module_lookup <- make_gene_lookup(panel_gene_sets)
heatmap_selected_genes <- unlist(heatmap_gene_sets, use.names = FALSE)
panel_selected_genes <- unlist(panel_gene_sets, use.names = FALSE)
selected_genes <- union(heatmap_selected_genes, panel_selected_genes)

gene_key <- tpm_table |>
  filter(gene_name %in% selected_genes) |>
  distinct(gene_name, .keep_all = TRUE)

missing_genes <- setdiff(selected_genes, gene_key$gene_name)
if (length(missing_genes) > 0) {
  stop("Missing genes in TPM table: ", paste(missing_genes, collapse = ", "))
}

gene_key <- gene_key |>
  mutate(gene_name = factor(gene_name, levels = selected_genes)) |>
  arrange(gene_name)

plot_time_levels <- sample_info |>
  distinct(time, time_hours) |>
  arrange(time_hours) |>
  pull(time)

expr_long <- gene_key |>
  select(gene_id, gene_name, all_of(sample_columns)) |>
  pivot_longer(
    cols = all_of(sample_columns),
    names_to = "sample",
    values_to = "tpm"
  ) |>
  mutate(log2_tpm = log2(tpm + 1)) |>
  left_join(sample_info, by = "sample")

plot_expr_long <- expr_long |>
  add_module_metadata(panel_module_lookup, names(module_titles), panel_selected_genes, pattern_labels)

heatmap_expr_long <- expr_long |>
  add_module_metadata(heatmap_module_lookup, names(module_titles), heatmap_selected_genes, pattern_labels)

heatmap_means <- summarise_mean_expr(
  heatmap_expr_long,
  c("module", "pattern", "gene_name", "infection", "age", "time", "time_hours"),
  value_col = "tpm"
) |>
  mutate(
    module = factor(module, levels = names(module_titles)),
    infection = factor(infection, levels = infection_levels),
    age = factor(age, levels = age_levels)
  )

heatmap_df_virus <- heatmap_means |>
  filter(module == "virus") |>
  build_ratio_heatmap(
    panel_var = "age",
    panel_levels = age_levels,
    ratio_var = "infection",
    numerator_level = "Flu",
    denominator_level = "COVID"
  )

heatmap_df_age <- heatmap_means |>
  filter(module == "age") |>
  build_ratio_heatmap(
    panel_var = "infection",
    panel_levels = c("Flu", "COVID"),
    ratio_var = "age",
    numerator_level = "Aged",
    denominator_level = "Young"
  )

p_heatmap_virus <- make_heatmap_plot(
  heatmap_df = heatmap_df_virus,
  module_levels = unlist(heatmap_gene_sets$virus, use.names = FALSE),
  module_title = module_titles[["virus"]],
  heatmap_palette = heatmap_palette,
  y_label_margin_right = 1,
  plot_title_size = plot_title_size,
  strip_text_size = strip_text_size,
  heatmap_x_text_size = heatmap_x_text_size,
  heatmap_y_text_size = heatmap_y_text_size,
  legend_text_size = legend_text_size
) +
  theme(legend.position = "none")

p_heatmap_age <- make_heatmap_plot(
  heatmap_df = heatmap_df_age,
  module_levels = unlist(heatmap_gene_sets$age, use.names = FALSE),
  module_title = module_titles[["age"]],
  heatmap_palette = heatmap_palette,
  y_label_margin_right = 3,
  plot_title_size = plot_title_size,
  strip_text_size = strip_text_size,
  heatmap_x_text_size = heatmap_x_text_size,
  heatmap_y_text_size = heatmap_y_text_size,
  legend_text_size = legend_text_size
) +
  theme(legend.position = "none")

p_heatmap_scale_legend <- make_heatmap_scale_legend_plot(
  heatmap_palette = heatmap_palette,
  legend_text_size = legend_text_size
)

p_virus <- make_module_plot(
  module_name = "virus",
  expr_long = plot_expr_long,
  module_levels = unlist(panel_gene_sets$virus, use.names = FALSE),
  module_title = module_titles[["virus"]],
  condition_palette = condition_palette,
  age_palette = age_palette,
  plot_time_levels = plot_time_levels,
  axis_title_size = axis_title_size,
  axis_text_size = axis_text_size,
  legend_text_size = legend_text_size,
  plot_title_size = plot_title_size,
  strip_text_size = strip_text_size
)

p_age <- make_module_plot(
  module_name = "age",
  expr_long = plot_expr_long,
  module_levels = unlist(panel_gene_sets$age, use.names = FALSE),
  module_title = module_titles[["age"]],
  condition_palette = condition_palette,
  age_palette = age_palette,
  plot_time_levels = plot_time_levels,
  axis_title_size = axis_title_size,
  axis_text_size = axis_text_size,
  legend_text_size = legend_text_size,
  plot_title_size = plot_title_size,
  strip_text_size = strip_text_size
)

top_heatmaps <- p_heatmap_virus | p_heatmap_age
heatmap_legends <- wrap_elements(
  full = plot_spacer() |
    p_heatmap_scale_legend |
    plot_spacer() +
    plot_layout(widths = c(0.8, 1, 0.8)),
  ignore_tag = TRUE
)
bottom_row <- p_virus | p_age

figure4 <- top_heatmaps / heatmap_legends / bottom_row +
  plot_layout(heights = c(1.02, 0.18, 1.3)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 14, face = "bold"))

save_figure_pdf(figure4, output_pdf, figure_width, figure_height)
