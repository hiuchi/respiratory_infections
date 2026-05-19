library(tidyverse)
library(patchwork)
library(scales)

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

resolve_named_palette <- function(levels, palette) {
  palette <- palette[intersect(names(palette), levels)]
  missing_levels <- setdiff(levels, names(palette))
  if (length(missing_levels) > 0) {
    extra_cols <- setNames(hue_pal()(length(missing_levels)), missing_levels)
    palette <- c(palette, extra_cols)
  }
  palette[levels]
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

project_dir <- "/Users/hiuchi/Dropbox/Research/covid/260309_analysis"
data_dir <- file.path(project_dir, "multiqc_data")
output_pdf <- build_figure_output(project_dir, "figure2.pdf")

condition_order <- c("Mock", "Flu", "COVID")
series_order <- c("Young", "Aged")
position_gaps <- list(time = 1, series = 2, condition = 5)

figure_width <- 8.27
figure_height <- 11.69

read_multiqc_tsv <- function(filename) {
  readr::read_tsv(
    file.path(data_dir, filename),
    na = c("", "NA", ".", "null"),
    show_col_types = FALSE,
    name_repair = "minimal"
  )
}

to_numeric <- function(x) {
  readr::parse_number(as.character(x))
}

preferred_levels <- function(values, preferred) {
  c(
    intersect(preferred, unique(values)),
    sort(setdiff(unique(values), preferred))
  )
}

assign_panel_positions <- function(df, time_gap = 1, series_gap = 2) {
  df <- arrange(df, series, time_hours, replicate_num, Sample)
  positions <- numeric(nrow(df))
  current_position <- 1

  for (i in seq_len(nrow(df))) {
    positions[i] <- current_position
    if (i < nrow(df)) {
      same_series <- df$series[i] == df$series[i + 1]
      same_time <- same_series && df$time_label[i] == df$time_label[i + 1]
      step <- if (same_time) {
        1
      } else if (same_series) {
        1 + time_gap
      } else {
        1 + series_gap
      }
      current_position <- current_position + step
    }
  }

  df$local_x_position <- positions
  df
}

build_axis_breaks <- function(sample_info) {
  sample_info |>
    distinct(condition, series, time_label, x_position) |>
    arrange(condition, series, time_label, x_position) |>
    group_by(condition, series, time_label) |>
    summarise(
      x_break = mean(x_position),
    .groups = "drop"
  ) |>
  mutate(axis_label = paste(series, time_label, sep = "\n"))
}

panel_x_expand <- expansion(mult = c(0.01, 0.01))

panel_x_scale_hidden <- scale_x_continuous(expand = panel_x_expand)

panel_x_scale_grouped <- function(axis_breaks) {
  scale_x_continuous(
    breaks = axis_breaks$x_break,
    labels = axis_breaks$axis_label,
    expand = panel_x_expand
  )
}

theme_hide_x <- theme(
  panel.grid.major.x = element_blank(),
  axis.text.x = element_blank(),
  axis.ticks.x = element_blank()
)

theme_grouped_x <- theme(
  panel.grid.major.x = element_blank(),
  axis.text.x = element_text(
    angle = 0,
    hjust = 0.5,
    vjust = 1,
    size = 8,
    lineheight = 0.9
  )
)

build_sample_info <- function(samples) {
  sample_info <- parse_sample_metadata(unique(samples)) |>
    mutate(
      Sample = sample_id,
      series_code = coalesce(series_code, "other"),
      condition = coalesce(condition, "Unknown"),
      time_label = coalesce(time_label, "NA"),
      replicate = coalesce(replicate, sample_id),
      series = coalesce(series, "Other"),
      time_hours = if_else(is.finite(time_hours), time_hours, Inf),
      replicate_num = coalesce(replicate_num, readr::parse_number(replicate))
    ) |>
    select(
      Sample,
      series_code,
      condition,
      time_label,
      replicate,
      series,
      time_hours,
      replicate_num
    )

  condition_levels <- preferred_levels(sample_info$condition, condition_order)
  series_levels <- preferred_levels(sample_info$series, series_order)

  sample_info <- sample_info |>
    mutate(
      condition = factor(condition, levels = unique(condition_levels)),
      series = factor(series, levels = unique(series_levels))
    ) |>
    arrange(condition, series, time_hours, replicate_num, Sample) |>
    group_by(condition) |>
    group_modify(
      ~ assign_panel_positions(
        .x,
        time_gap = position_gaps$time,
        series_gap = position_gaps$series
      )
    ) |>
    ungroup()

  condition_offsets <- sample_info |>
    group_by(condition) |>
    summarise(local_max = max(local_x_position), .groups = "drop") |>
    mutate(offset = lag(cumsum(local_max + position_gaps$condition), default = 0))

  sample_info <- sample_info |>
    left_join(condition_offsets, by = "condition") |>
    mutate(x_position = local_x_position + offset)

  sample_order <- sample_info |>
    arrange(condition, series, time_hours, replicate_num, Sample) |>
    pull(Sample)

  sample_info |>
    mutate(Sample = factor(Sample, levels = sample_order))
}

join_sample_info <- function(df, sample_info) {
  df |>
    left_join(sample_info, by = "Sample") |>
    mutate(Sample = factor(Sample, levels = levels(sample_info$Sample)))
}

condition_palette <- figure_palettes$infection
origin_palette <- c(
  Exonic = "#1B9E77",
  Intronic = "#D95F02",
  Intergenic = "#7570B3"
)

base_theme <- make_figure_theme(
  base_size = 12,
  plot_title_size = 13,
  axis_title_size = 12,
  axis_text_size = 11,
  legend_text_size = 11,
  strip_text_size = 12
) +
  theme(
    legend.position = "top",
    plot.subtitle = element_text(size = 11)
  )

theme_set(base_theme)

star_summary <- read_multiqc_tsv("star_summary_table.tsv") |>
  transmute(
    Sample,
    total_reads_m = to_numeric(`Total reads`),
    uniquely_mapped_percent = to_numeric(`Uniq aligned`)
  ) |>
  filter(!is.na(total_reads_m), !is.na(uniquely_mapped_percent))

genomic_origin <- read_multiqc_tsv("qualimap_genomic_origin.tsv") |>
  transmute(
    Sample,
    Exonic = to_numeric(Exonic),
    Intronic = to_numeric(Intronic),
    Intergenic = to_numeric(Intergenic)
  ) |>
  filter(!if_any(-Sample, is.na))

general_stats_raw <- read_multiqc_tsv("general_stats_table.tsv")
general_stats <- tibble(
  Sample = general_stats_raw[["Sample"]],
  dup_int = to_numeric(general_stats_raw[["dupInt"]]),
  duplication_percent = to_numeric(general_stats_raw[["Duplication"]])
) |>
  filter(!str_detect(Sample, " Read [12]$")) |>
  filter(!if_all(c(dup_int, duplication_percent), is.na))

coverage_profile <- read_multiqc_tsv("qualimap_gene_coverage_profile.tsv") |>
  pivot_longer(
    cols = -Sample,
    names_to = "gene_body_percent",
    values_to = "coverage"
  ) |>
  mutate(
    position = to_numeric(gene_body_percent),
    coverage = to_numeric(coverage)
  ) |>
  filter(!is.na(position), !is.na(coverage)) |>
  group_by(Sample) |>
  mutate(relative_coverage = coverage / mean(coverage, na.rm = TRUE)) |>
  ungroup()

sample_info <- build_sample_info(
  c(
    star_summary$Sample,
    genomic_origin$Sample,
    general_stats$Sample,
    coverage_profile$Sample
  )
)

star_summary <- join_sample_info(star_summary, sample_info)
genomic_origin <- join_sample_info(genomic_origin, sample_info)
general_stats <- join_sample_info(general_stats, sample_info)
coverage_profile <- join_sample_info(coverage_profile, sample_info)

condition_palette <- resolve_named_palette(levels(sample_info$condition), condition_palette)

origin_long <- genomic_origin |>
  pivot_longer(
    cols = c(Exonic, Intronic, Intergenic),
    names_to = "origin",
    values_to = "percent"
  ) |>
  mutate(origin = factor(origin, levels = c("Exonic", "Intronic", "Intergenic")))

origin_axis_breaks <- build_axis_breaks(sample_info)

coverage_summary <- coverage_profile |>
  group_by(condition, series, position) |>
  summarise(mean_relative_coverage = mean(relative_coverage, na.rm = TRUE), .groups = "drop")

complexity_correlation <- cor(
  general_stats$duplication_percent,
  general_stats$dup_int,
  use = "complete.obs"
)

p_depth <- ggplot(star_summary, aes(x_position, total_reads_m, fill = condition)) +
  geom_col(width = 0.85, colour = NA) +
  facet_grid(cols = vars(condition), scales = "free_x", space = "free_x") +
  scale_fill_manual(values = condition_palette) +
  panel_x_scale_hidden +
  guides(fill = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.03))) +
  labs(
    title = "A Sequencing depth",
    y = "Total reads (millions)",
    x = NULL
  ) +
  theme(legend.position = "none") +
  theme_hide_x

p_alignment <- ggplot(star_summary, aes(x_position, uniquely_mapped_percent, colour = condition)) +
  geom_point(size = 1.8, alpha = 0.9) +
  facet_grid(cols = vars(condition), scales = "free_x", space = "free_x") +
  scale_colour_manual(values = condition_palette) +
  panel_x_scale_hidden +
  guides(colour = "none") +
  scale_y_continuous(
    limits = c(0, 100),
    labels = label_number(suffix = "%")
  ) +
  labs(
    title = "B Alignment quality",
    y = "Uniquely mapped reads (%)",
    x = NULL
  ) +
  theme(legend.position = "none") +
  theme_hide_x

p_origin <- ggplot(origin_long, aes(x_position, percent, fill = origin)) +
  geom_col(width = 0.9, colour = NA) +
  facet_grid(cols = vars(condition), scales = "free_x", space = "free_x") +
  scale_fill_manual(values = origin_palette) +
  panel_x_scale_grouped(origin_axis_breaks) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.02)),
    labels = label_number(suffix = "%")
  ) +
  labs(
    title = "C Genomic origin of reads",
    y = "Fraction of aligned reads",
    x = NULL
  ) +
  theme_grouped_x

p_complexity <- ggplot(general_stats, aes(duplication_percent, dup_int, colour = condition, shape = series)) +
  geom_smooth(
    data = general_stats,
    inherit.aes = FALSE,
    aes(duplication_percent, dup_int, group = 1),
    method = "lm",
    formula = y ~ x,
    se = FALSE,
    linewidth = 0.5,
    linetype = 2,
    colour = "#6B6B6B"
  ) +
  geom_point(size = 2.3, alpha = 0.9) +
  scale_colour_manual(values = condition_palette) +
  scale_x_continuous(labels = label_number(suffix = "%")) +
  labs(
    title = "D Library complexity",
    subtitle = paste0("Duplication vs dupInt, Pearson r = ", sprintf("%.2f", complexity_correlation)),
    x = "Duplication (%)",
    y = "dupInt"
  )

p_coverage <- ggplot() +
  geom_hline(yintercept = 1, linewidth = 0.4, linetype = 2, colour = "#A0A0A0") +
  geom_line(
    data = coverage_profile,
    aes(position, relative_coverage, group = Sample, colour = condition),
    linewidth = 0.25,
    alpha = 0.12
  ) +
  geom_line(
    data = coverage_summary,
    aes(position, mean_relative_coverage, colour = condition, linetype = series),
    linewidth = 1
  ) +
  scale_colour_manual(values = condition_palette) +
  guides(colour = "none", linetype = "none") +
  scale_x_continuous(
    breaks = c(0, 25, 50, 75, 99),
    labels = c("0", "25", "50", "75", "100")
  ) +
  labs(
    title = "E Gene coverage profile",
    subtitle = "Qualimap profile normalized to the per-sample mean coverage",
    x = "Gene body position (%)",
    y = "Relative coverage"
  ) +
  theme(legend.position = "none")

figure2 <- (p_depth / p_alignment / p_origin / (p_complexity | p_coverage)) +
  plot_layout(heights = c(1.1, 1.1, 1.4, 1.1), guides = "collect") &
  theme(legend.position = "top")

save_figure_pdf(figure2, output_pdf, figure_width, figure_height)
