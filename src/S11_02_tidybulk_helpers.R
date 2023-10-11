#!/usr/bin/env Rscript

tidybulk_helpers <- new.env()

tidybulk_helpers$plot_normalization <- function(tidybulk_object) {
  tidybulk_object %>%
    tidyr::pivot_longer(
      c(counts, counts_scaled, counts_scaled_adjusted),
      values_to = "count",
      names_to = "normalization"
    ) %>%
    ggplot2::ggplot(aes(count + 1, color = batch_effect)) +
    ggplot2::geom_density(na.rm = TRUE) +
    ggplot2::scale_x_log10() +
    ggplot2::facet_grid(~normalization)
}


tidybulk_helpers$ggsave_a4 <- function(plot, filename) {
  ggplot2::ggsave(
    plot = plot, filename = filename, height = 210, width = 297, units = "mm"
  )
}


tidybulk_helpers$ggpairs_pca <- function(data, column) {
  column_q <- substitute(column)

  GGally::ggpairs(
    data = data,
    legend = c(2, 1),
    title = str_glue("PCA by {column_q}"),
    columns = 1:10,
    ggplot2::aes(
      color = !!column_q,
      sequencing_batch = sequencing_batch, lab_batch = lab_batch, ribo_batch = ribo_batch,
      sample = sample, breed = breed, sex = sex, day = day,
      treatment = treatment, campylo = has_campylobacter
    ),
    upper = NULL, diag = NULL
  )
}


tidybulk_helpers$annotate_differential_expression <- function(tidybulk_de) {
  tidybulk_de %>%
    tidybulk::keep_abundant() %>%
    tidybulk::pivot_transcript() %>%
    dplyr::mutate(
      significant = FDR < 0.05 & abs(logFC) > 1.2,
      direction = sign(logFC),
      label = dplyr::if_else(significant, gene, NA_character_),
    )
}

tidybulk_helpers$plot_ma <- function(tidybulk_de_annotated) {
  tidybulk_de_annotated %>%
    ggplot2::ggplot(
      ggplot2::aes(x = logCPM, y = logFC, label = label)
    ) +
    ggplot2::geom_point(
      ggplot2::aes(color = significant, size = significant, alpha = significant),
      na.rm = TRUE
    ) +
    ggplot2::scale_color_manual(values = c("black", "red")) +
    ggplot2::scale_size_manual(values = c(1, 2)) +
    ggplot2::scale_alpha_manual(values = c(0.1, 1)) +
    ggrepel::geom_text_repel(na.rm = TRUE)
}

tidybulk_helpers$plot_volcano <- function(tidybulk_de_annotated) {
  tidybulk_de_annotated %>%
    ggplot2::ggplot(aes(x = logFC, y = FDR, label = label)) +
    ggplot2::geom_point(
      ggplot2::aes(color = significant, size = significant, alpha = significant),
      na.rm = TRUE
    ) +
    ggplot2::scale_color_manual(values = c("black", "red")) +
    ggplot2::scale_size_manual(values = c(1, 2)) +
    ggplot2::scale_alpha_manual(values = c(0.1, 1)) +
    ggplot2::scale_y_continuous(trans = "log10_reverse", na.value = NA) +
    ggrepel::geom_text_repel(na.rm = TRUE)
}

attach(tidybulk_helpers, name = "tidybulk_helpers")
rm(tidybulk_helpers)
