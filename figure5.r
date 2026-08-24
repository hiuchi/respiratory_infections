library(DESeq2)
library(patchwork)
library(tidyverse)

figure_strip_fill <- "#F4F4F4"
figure_strip_border <- "#D0D0D0"
a4_landscape_width <- 11.69
a4_landscape_height <- 8.27

build_output_dir <- function(project_dir, subdir) {
  output_dir <- file.path(project_dir, "output", subdir)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_dir
}

build_data_output <- function(project_dir, filename) {
  file.path(build_output_dir(project_dir, "data"), filename)
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

make_figure_theme <- function(
  base_size = 11,
  plot_title_size = base_size + 1,
  axis_title_size = base_size,
  axis_text_size = base_size - 1,
  legend_text_size = base_size - 1,
  strip_text_size = base_size,
  plot_title_position = "plot"
) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      legend.title = element_blank(),
      strip.background = element_rect(fill = figure_strip_fill, colour = figure_strip_border),
      strip.text = element_text(size = strip_text_size, face = "bold"),
      axis.title = element_text(size = axis_title_size),
      axis.text = element_text(size = axis_text_size),
      legend.text = element_text(size = legend_text_size),
      plot.title = element_text(size = plot_title_size, face = "bold"),
      plot.title.position = plot_title_position
    )
}

extract_gff3_attribute <- function(attributes, key) {
  str_match(attributes, paste0("(?:^|;)", key, "=([^;]+)"))[, 2]
}

read_gene_annotation <- function(gff3_path) {
  read_tsv(
    gff3_path,
    col_names = c(
      "seqname", "source", "feature", "start", "end",
      "score", "strand", "phase", "attributes"
    ),
    col_types = cols(.default = col_character()),
    comment = "#",
    show_col_types = FALSE
  ) |>
    filter(feature == "gene") |>
    transmute(
      gene_id = extract_gff3_attribute(attributes, "gene_id"),
      ensembl_id = str_remove(gene_id, "\\..*$"),
      gene_name = extract_gff3_attribute(attributes, "gene_name"),
      gene_type = coalesce(
        extract_gff3_attribute(attributes, "gene_type"),
        extract_gff3_attribute(attributes, "gene_biotype")
      )
    ) |>
    filter(!is.na(ensembl_id)) |>
    distinct(ensembl_id, .keep_all = TRUE)
}

parse_sample_info <- function(samples) {
  parts <- str_split_fixed(samples, "-", 4)

  tibble(
    sample = samples,
    age = case_when(
      parts[, 1] == "y" ~ "Young",
      parts[, 1] == "a" ~ "Aged",
      TRUE ~ NA_character_
    ),
    infection = parts[, 2],
    day = parts[, 3],
    replicate = parts[, 4]
  ) |>
    filter(!is.na(age), infection != "", day != "") |>
    mutate(
      age = factor(age, levels = c("Young", "Aged")),
      infection = factor(
        infection,
        levels = c("Mock", "Flu", "COVID"),
        labels = c("Mock", "Influenza", "COVID-19")
      ),
      day = factor(day, levels = c("6h", "1d", "2d", "4d", "6d", "8d", "10d"))
    )
}

make_group_result <- function(normalized_counts, sample_info, annotation, age_value, infection_value, day_value) {
  target_samples <- sample_info |>
    filter(age == age_value, infection == infection_value, day == day_value) |>
    pull(sample)
  mock_samples <- sample_info |>
    filter(age == age_value, infection == "Mock", day == "6h") |>
    pull(sample)

  if (length(target_samples) == 0 || length(mock_samples) == 0) {
    return(tibble())
  }

  target_mean <- rowMeans(normalized_counts[, target_samples, drop = FALSE])
  mock_mean <- rowMeans(normalized_counts[, mock_samples, drop = FALSE])

  tibble(
    gene_id = rownames(normalized_counts),
    ensembl_id = str_remove(gene_id, "\\..*$"),
    age = age_value,
    infection = infection_value,
    day = day_value,
    mock_mean = mock_mean,
    infection_mean = target_mean,
    log2FoldChange = log2(infection_mean / mock_mean)
  ) |>
    left_join(annotation, by = c("gene_id", "ensembl_id")) |>
    mutate(
      type = case_when(
        gene_type == "protein_coding" ~ "mRNA",
        gene_type == "lncRNA" ~ "lncRNA",
        TRUE ~ "Others"
      )
    ) |>
    filter(mock_mean >= 1, is.finite(log2FoldChange))
}

format_p_value <- function(p_value) {
  case_when(
    p_value < 0.0005 ~ "***",
    p_value < 0.005 ~ "**",
    p_value < 0.05 ~ "*",
    TRUE ~ "N.S."
  )
}

build_p_value_labels <- function(plot_data) {
  plot_data |>
    group_by(age, day) |>
    group_modify(function(data, key) {
      group_a <- data$log2FoldChange[data$type == "mRNA"]
      group_b <- data$log2FoldChange[data$type == "lncRNA"]
      p_value <- wilcox.test(group_a, group_b, exact = FALSE)$p.value

      tibble(
        x = 1.5,
        y = 2.55,
        label = format_p_value(p_value)
      )
    }) |>
    ungroup()
}

plot_timecourse <- function(plot_data, panel_label) {
  p_value_labels <- build_p_value_labels(plot_data)

  ggplot(plot_data, aes(x = type, y = log2FoldChange, colour = type)) +
    geom_boxplot(outlier.shape = NA, width = 0.72) +
    geom_text(
      data = p_value_labels,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      size = 4,
      fontface = "bold"
    ) +
    facet_grid(age ~ day, drop = FALSE) +
    coord_cartesian(ylim = c(-2.75, 2.75)) +
    scale_colour_manual(values = c(mRNA = "#4C78A8", lncRNA = "#D65F5F")) +
    labs(x = NULL, y = "Log2 fold change", title = panel_label) +
    make_figure_theme(
      base_size = 12,
      plot_title_size = 18,
      axis_title_size = 12,
      axis_text_size = 10,
      legend_text_size = 11,
      strip_text_size = 12,
      plot_title_position = "panel"
    ) +
    theme(
      axis.text.x = element_blank(),
      legend.position = "bottom",
      plot.title = element_text(size = 18, face = "bold", margin = margin(b = 4))
    )
}

project_dir <- "/path/to/data"
input_dir <- file.path(project_dir, "res", "star_salmon", "deseq2_qc")
annotation_gff3 <- file.path(project_dir, "files", "mouse_virus.gff3")
plot_output_dir <- build_output_dir(project_dir, "plots")

load(file.path(input_dir, "deseq2.dds.RData"))

sample_info <- parse_sample_info(colnames(dds))
normalized_counts <- counts(dds, normalized = TRUE)
annotation <- read_gene_annotation(annotation_gff3)

group_specs <- sample_info |>
  filter(infection %in% c("Influenza", "COVID-19")) |>
  distinct(age, infection, day) |>
  arrange(infection, age, day) |>
  transmute(
    age_value = as.character(age),
    infection_value = as.character(infection),
    day_value = as.character(day)
  )

results <- purrr::pmap_dfr(
  group_specs,
  make_group_result,
  normalized_counts = normalized_counts,
  sample_info = sample_info,
  annotation = annotation
) |>
  mutate(
    age = factor(age, levels = c("Young", "Aged")),
    infection = factor(infection, levels = c("Mock", "Influenza", "COVID-19")),
    day = factor(day, levels = c("6h", "1d", "2d", "4d", "6d", "8d", "10d")),
    type = factor(type, levels = c("mRNA", "lncRNA", "Others"))
  )

saveRDS(results, build_data_output(project_dir, "lncRNA_timecourse_results.rds"))
write_csv(results, build_data_output(project_dir, "lncRNA_timecourse_results.csv"))

plot_data <- results |>
  filter(type %in% c("mRNA", "lncRNA"))

flu_plot <- plot_timecourse(plot_data |> filter(infection == "Influenza", !is.na(day)), "A  Influenza")
covid_plot <- plot_timecourse(plot_data |> filter(infection == "COVID-19", !is.na(day)), "B  COVID-19")

figure5 <- ((flu_plot / covid_plot) + plot_layout(guides = "collect")) &
  theme(legend.position = "bottom")

save_figure_pdf(
  figure5,
  file.path(plot_output_dir, "figure5.pdf"),
  figure_width = a4_landscape_width,
  figure_height = a4_landscape_height
)
