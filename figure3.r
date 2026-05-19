library(DESeq2)
library(cowplot)
library(ggdendro)
library(matrixStats)
library(patchwork)
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

input_dir <- "/Users/hiuchi/Dropbox/Research/covid/260305_analysis/res/star_salmon/deseq2_qc"
project_dir <- "/Users/hiuchi/Dropbox/Research/covid/260309_analysis"
output_pdf <- build_figure_output(project_dir, "figure3.pdf")
figure_width <- 15
figure_height <- 11.2

infection_palette <- figure_palettes$infection
age_palette <- figure_palettes$age
pca_shape_values <- c(Young = 16, Aged = 17)
distance_palette <- c(
  "#B2182B",
  "#D6604D",
  "#F4A582",
  "#FDDBC7",
  "#D1E5F0",
  "#92C5DE",
  "#4393C3",
  "#2166AC"
)
legend_key_height <- grid::unit(9, "pt")
heatmap_label_size <- 6

theme_set(
  make_figure_theme(
    base_size = 11,
    plot_title_size = 11,
    axis_title_size = 11,
    axis_text_size = 10,
    legend_text_size = 10,
    strip_text_size = 11
  )
)

parse_sample_info <- function(samples) {
  parse_sample_metadata(samples) |>
    filter(
      is_valid,
      series %in% c("Young", "Aged"),
      condition %in% c("Mock", "Flu", "COVID"),
      !is.na(replicate_num)
    ) |>
    transmute(
      sample = sample_id,
      age = factor(series, levels = c("Young", "Aged")),
      infection = factor(condition, levels = c("Mock", "Flu", "COVID")),
      time = time_label,
      replicate = as.integer(replicate_num),
      time_hours
    )
}

make_sample_labels <- function(sample_info) {
  sample_info |>
    transmute(
      sample,
      label = sprintf(
        "%s-%s",
        recode(as.character(age), Young = "Y", Aged = "A"),
        time
      )
    )
}

make_label_lookup <- function(label_df) {
  setNames(label_df$label, label_df$sample)
}

make_pca_axis_labels <- function(variance_explained) {
  c(
    x = sprintf("PC1 (%.1f%%)", 100 * variance_explained[1]),
    y = sprintf("PC2 (%.1f%%)", 100 * variance_explained[2])
  )
}

make_square_pca_theme <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      aspect.ratio = 1,
      legend.position = "none",
      plot.title = element_blank(),
      plot.subtitle = element_blank()
    )
}

wrap_patchwork_panel <- function(plot) {
  panel_grob <- local({
    grDevices::pdf(NULL)
    on.exit(invisible(grDevices::dev.off()))

    patchworkGrob(plot)
  })

  wrap_elements(
    full = panel_grob,
    clip = FALSE
  )
}

make_bottom_legend_theme <- function(text_size, title_size = text_size, box = NULL) {
  theme_args <- list(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.text = element_text(size = text_size),
    legend.title = element_text(size = title_size),
    legend.key.height = legend_key_height,
    legend.margin = margin(0, 0, 0, 0)
  )

  if (!is.null(box)) {
    theme_args$legend.box <- box
  }

  do.call(theme, theme_args)
}

make_annotation_df <- function(sample_info) {
  annotation_df <- sample_info |>
    transmute(
      sample,
      Infection = infection,
      Age = age
    ) |>
    as.data.frame()
  rownames(annotation_df) <- annotation_df$sample
  annotation_df$sample <- NULL
  annotation_df
}

make_dend_bottom_df <- function(dend_meta) {
  bind_rows(
    dend_meta |>
      transmute(
        x,
        y = 2,
        fill = unname(age_palette[as.character(age)])
      ),
    dend_meta |>
      transmute(
        x,
        y = 1,
        fill = unname(infection_palette[as.character(infection)])
      )
  )
}

make_correlation_panel <- function(heatmap_gtable) {
  grDevices::pdf(NULL)
  on.exit(invisible(grDevices::dev.off()))

  ggdraw() +
    draw_grob(heatmap_gtable, x = 0, y = 0, width = 1, height = 1) +
    draw_label("Age", x = 0.022, y = 0.065, hjust = 0, size = 8) +
    draw_label("Infection", x = 0.022, y = 0.03, hjust = 0, size = 8)
}

load(file.path(input_dir, "deseq2.dds.RData"))

sample_info <- parse_sample_info(colnames(dds))
sample_labels <- make_sample_labels(sample_info)
sample_label_lookup <- make_label_lookup(sample_labels)
annotation_colors <- list(
  Infection = infection_palette,
  Age = age_palette
)

vsd <- vst(dds, blind = TRUE)
vst_mat <- assay(vsd)

ntop <- min(500L, nrow(vst_mat))
top_idx <- order(rowVars(vst_mat), decreasing = TRUE)[seq_len(ntop)]
pca_fit <- prcomp(t(vst_mat[top_idx, ]), center = TRUE, scale. = FALSE)
variance_explained <- (pca_fit$sdev^2) / sum(pca_fit$sdev^2)
pca_axis_labels <- make_pca_axis_labels(variance_explained)

pca_scores <- as_tibble(pca_fit$x[, 1:3], rownames = "sample") |>
  left_join(sample_info, by = "sample")

time_levels <- sample_info |>
  distinct(time, time_hours) |>
  arrange(time_hours)

time_palette <- setNames(
  grDevices::hcl.colors(nrow(time_levels), palette = "Temps"),
  time_levels$time
)

cor_mat <- cor(vst_mat, method = "pearson")
dist_mat <- 1 - cor_mat
distance_limits <- c(0, ceiling(max(dist_mat, na.rm = TRUE) * 100) / 100)
hc <- hclust(as.dist(dist_mat), method = "average")

dend_data <- local({
  grDevices::pdf(NULL)
  on.exit(invisible(grDevices::dev.off()))

  ggdendro::dendro_data(as.dendrogram(hc), type = "rectangle")
})
dend_segments <- segment(dend_data)
dend_labels <- as_tibble(dend_data$labels) |>
  transmute(sample = label, x = x)

dend_meta <- sample_info |>
  inner_join(dend_labels, by = "sample") |>
  mutate(label = unname(sample_label_lookup[sample])) |>
  arrange(x)

dend_x_limits <- c(0.5, max(dend_labels$x) + 0.5)

p_pca12 <- ggplot(pca_scores, aes(PC1, PC2, colour = infection, shape = age)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "#C8C8C8") +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "#C8C8C8") +
  geom_point(size = 2.2, alpha = 0.9) +
  scale_colour_manual(values = infection_palette, name = NULL) +
  scale_shape_manual(values = pca_shape_values, name = NULL) +
  labs(
    x = pca_axis_labels["x"],
    y = pca_axis_labels["y"]
  ) +
  make_square_pca_theme() +
  theme(plot.margin = margin(6, 6, 6, 6))

pca_scores_time <- pca_scores |>
  mutate(time = factor(time, levels = time_levels$time))

p_pca_time <- ggplot(
  pca_scores_time,
  aes(PC1, PC2, colour = time, shape = age)
) +
  geom_point(size = 2.2, alpha = 0.9) +
  scale_colour_manual(
    values = time_palette,
    breaks = time_levels$time,
    labels = time_levels$time,
    name = NULL
  ) +
  scale_shape_manual(values = pca_shape_values, name = NULL) +
  labs(
    x = pca_axis_labels["x"],
    y = pca_axis_labels["y"]
  ) +
  make_square_pca_theme(base_size = 11) +
  theme(
    plot.margin = margin(6, 6, 6, 6)
  )

p_pca12_legend <- p_pca12 +
  guides(
    colour = guide_legend(order = 1, title = NULL, nrow = 1, override.aes = list(shape = 16, size = 3)),
    shape = guide_legend(order = 2, title = NULL, nrow = 1)
  ) +
  make_bottom_legend_theme(text_size = 8, box = "vertical")

p_pca_time_legend <- p_pca_time +
  guides(
    colour = guide_legend(
      order = 1,
      title = NULL,
      nrow = 2,
      byrow = TRUE,
      override.aes = list(shape = 16, size = 2.8)
    ),
    shape = guide_legend(order = 2, title = NULL, nrow = 1)
  ) +
  make_bottom_legend_theme(text_size = 7, title_size = 8, box = "vertical")

a_panel_raw <- (p_pca12_legend | p_pca_time_legend) + plot_layout(widths = c(1, 1))

a_panel <- wrap_patchwork_panel(a_panel_raw)

annotation_df <- make_annotation_df(sample_info)

dist_breaks <- seq(distance_limits[1], distance_limits[2], length.out = 101)
cor_heatmap <- local({
  grDevices::pdf(NULL)
  on.exit(invisible(grDevices::dev.off()))

  pheatmap::pheatmap(
    mat = dist_mat,
    color = grDevices::colorRampPalette(distance_palette)(100),
    breaks = dist_breaks,
    cluster_rows = hc,
    cluster_cols = hc,
    annotation_row = annotation_df,
    annotation_col = annotation_df,
    annotation_colors = annotation_colors,
    show_rownames = TRUE,
    show_colnames = TRUE,
    labels_row = unname(sample_label_lookup[rownames(dist_mat)]),
    labels_col = unname(sample_label_lookup[colnames(dist_mat)]),
    fontsize = 8,
    fontsize_row = heatmap_label_size,
    fontsize_col = heatmap_label_size,
    angle_col = 90,
    border_color = "#D9D9D9",
    treeheight_row = 42,
    treeheight_col = 42,
    cellwidth = 4.1,
    cellheight = 4.5,
    annotation_names_row = FALSE,
    annotation_names_col = TRUE,
    silent = TRUE
  )
})
cor_heatmap$gtable$heights[5] <- grid::unit(78, "bigpts")

cor_panel <- make_correlation_panel(cor_heatmap$gtable)

dend_max_height <- max(c(dend_segments$y, dend_segments$yend))

p_dendro <- ggplot(dend_segments) +
  geom_segment(
    aes(x = x, y = y, xend = xend, yend = yend),
    linewidth = 0.35,
    colour = "#333333",
    lineend = "square"
  ) +
  scale_x_continuous(
    limits = dend_x_limits,
    expand = c(0, 0)
  ) +
  coord_cartesian(
    ylim = c(0, dend_max_height * 1.03),
    clip = "off"
  ) +
  labs(
    x = NULL,
    y = "Height"
  ) +
  theme_classic(base_size = 10) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(2, 4, 0, 4),
    plot.title = element_blank(),
    plot.subtitle = element_blank()
  )

dend_bottom_df <- make_dend_bottom_df(dend_meta)

p_dend_bottom <- ggplot() +
  geom_tile(
    data = dend_bottom_df,
    aes(x = x, y = y, fill = fill),
    width = 0.95,
    height = 0.85
  ) +
  scale_fill_identity() +
  scale_x_continuous(
    limits = dend_x_limits,
    breaks = dend_meta$x,
    labels = dend_meta$label,
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(0.4, 2.6),
    breaks = c(1, 2),
    labels = c("Infection", "Age"),
    expand = c(0, 0)
  ) +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 8.5) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6.1),
    axis.text.y = element_text(size = 7),
    axis.ticks.x = element_blank(),
    plot.margin = margin(0, 4, 2, 4)
  )

dend_panel_raw <- p_dendro / p_dend_bottom + plot_layout(heights = c(1, 0.42))
dend_panel <- wrap_patchwork_panel(dend_panel_raw)

top_row <- (a_panel | cor_panel) + plot_layout(widths = c(1.45, 2.05))

figure3 <- (top_row / dend_panel) +
  plot_layout(heights = c(1.12, 0.72)) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(size = 14, face = "bold")
  )

save_figure_pdf(figure3, output_pdf, figure_width, figure_height)
