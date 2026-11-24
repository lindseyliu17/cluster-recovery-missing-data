# =============================================================================
# 07_results.R
#
# PURPOSE: Aggregate simulation results and produce the tables and figures
#          for the final report.
#
# OUTPUT FILES
# ------------
# Tables (CSV, in output/tables/):
#   table1_ari_by_proximity_rate.csv   — mean ARI per proximity × miss rate
#   table2_ari_by_clustering_method.csv — mean ARI per proximity × method
#   table3_zero_overlap_pairs.csv       — pairwise no-overlap diagnostic
#
# Figures (PNG, in output/figures/):
#   fig1_continuous_ari_boxplot_rate.png, fig1_mixed_ari_boxplot_rate.png
#   fig2_continuous_ari_lineplot_rate.png, fig2_mixed_ari_lineplot_rate.png
#   fig3_continuous_ari_heatmap.png, fig3_mixed_ari_heatmap.png
#   fig4_continuous_silhouette_boxplot.png, fig4_mixed_silhouette_boxplot.png
# =============================================================================

# Readable display labels for the five proximity measures
PROXIMITY_LABELS <- c(
  euclidean_imputed  = "Eucl. (Imputed)",
  euclidean_pairwise = "Eucl. (Pairwise)",
  manhattan_imputed  = "Manh. (Imputed)",
  manhattan_pairwise = "Manh. (Pairwise)",
  gower              = "Gower"
)

# Preferred ordering for plots (Gower last, paired measures together)
PROXIMITY_ORDER <- c(
  "Eucl. (Imputed)", "Eucl. (Pairwise)",
  "Manh. (Imputed)", "Manh. (Pairwise)",
  "Gower"
)

DATA_TYPE_LABELS <- c(
  continuous = "Continuous",
  mixed      = "Mixed"
)

DATA_TYPE_FILE_STEMS <- c(
  Continuous = "continuous",
  Mixed      = "mixed"
)

FIGURE_THEME <- function(base_size = 12) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      plot.title      = ggplot2::element_text(face = "bold", hjust = 0.5),
      strip.background = ggplot2::element_rect(fill = "grey90", color = "grey35"),
      strip.text      = ggplot2::element_text(face = "bold"),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.box      = "vertical"
    )
}


#' Master function: build all tables and figures from the simulation results
#'
#' @param results  Data frame returned by run_full_simulation()
create_all_outputs <- function(results) {

  if (!"data_type" %in% names(results)) {
    results$data_type <- "continuous"
  }

  # ---- Prepare clean display labels -----------------------------------------
  results <- results %>%
    dplyr::mutate(
      miss_pct   = factor(paste0(miss_rate * 100, "%"),
                          levels = c("0%", "10%", "30%", "50%")),
      data_label = factor(DATA_TYPE_LABELS[data_type],
                          levels = c("Continuous", "Mixed")),
      prox_label = factor(PROXIMITY_LABELS[proximity],
                          levels = PROXIMITY_ORDER),
      clust_label = dplyr::recode(clustering,
        pam     = "PAM",
        average = "UPGMA"
      )
    )

  # ---- Tables ---------------------------------------------------------------
  cat("  Building tables...\n")
  make_table1_ari_by_rate(results)
  make_table2_ari_by_method(results)
  make_table3_zero_overlap_pairs(results)

  # ---- Figures --------------------------------------------------------------
  cat("  Building figures...\n")
  plot_fig1_ari_boxplot_rate(results)
  plot_fig2_ari_lineplot_rate(results)
  plot_fig3_ari_heatmap(results)
  plot_fig4_silhouette_boxplot(results)

  cat("  All outputs saved.\n")
}


# =============================================================================
# TABLE 3 — Mean proportion of observation pairs with no jointly observed
#           variable. This is defined for pairwise-complete Euclidean and
#           Manhattan distances, and for Gower. It is NA for imputed methods.
# =============================================================================
make_table3_zero_overlap_pairs <- function(results) {

  if (!"zero_pair_prop" %in% names(results)) {
    cat("    Skipped: zero_pair_prop not found in results.\n")
    return(invisible(NULL))
  }

  tbl <- results %>%
    dplyr::filter(!is.na(zero_pair_prop)) %>%
    dplyr::distinct(
      data_label, separation, mechanism, miss_pct, prox_label, replicate,
      zero_pair_prop, zero_pairs
    ) %>%
    dplyr::group_by(data_label, separation, mechanism, miss_pct, prox_label) %>%
    dplyr::summarise(
      mean_zero_pair_prop = round(mean(zero_pair_prop, na.rm = TRUE), 4),
      sd_zero_pair_prop   = round(sd(zero_pair_prop,   na.rm = TRUE), 4),
      mean_zero_pairs     = round(mean(zero_pairs,     na.rm = TRUE), 1),
      .groups             = "drop"
    ) %>%
    dplyr::rename(
      "Data Type" = data_label,
      "Separation" = separation,
      "Mechanism" = mechanism,
      "Missingness Rate" = miss_pct,
      "Proximity Measure" = prox_label,
      "Mean Zero-Overlap Pair Proportion" = mean_zero_pair_prop,
      "SD Zero-Overlap Pair Proportion" = sd_zero_pair_prop,
      "Mean Zero-Overlap Pair Count" = mean_zero_pairs
    )

  write.csv(tbl, "output/tables/table3_zero_overlap_pairs.csv",
            row.names = FALSE)
  cat("    Saved: output/tables/table3_zero_overlap_pairs.csv\n")
}


# =============================================================================
# TABLE 1 — Mean (SD) ARI by data type × proximity measure × missingness rate
#           (averaged over both clustering methods and both mechanisms)
# =============================================================================
make_table1_ari_by_rate <- function(results) {

  tbl <- results %>%
    dplyr::group_by(data_label, prox_label, miss_pct) %>%
    dplyr::summarise(
      mean_ari = round(mean(ari, na.rm = TRUE), 3),
      sd_ari   = round(sd(ari,   na.rm = TRUE), 3),
      .groups  = "drop"
    ) %>%
    dplyr::mutate(cell = paste0(mean_ari, " (", sd_ari, ")")) %>%
    dplyr::select(data_label, prox_label, miss_pct, cell) %>%
    tidyr::pivot_wider(names_from = miss_pct, values_from = cell) %>%
    dplyr::rename(
      "Data Type" = data_label,
      "Proximity Measure" = prox_label
    )

  write.csv(tbl, "output/tables/table1_ari_by_proximity_rate.csv",
            row.names = FALSE)
  cat("    Saved: output/tables/table1_ari_by_proximity_rate.csv\n")
}


# =============================================================================
# TABLE 2 — Mean ARI by data type × proximity measure × clustering method
#           (averaged over all miss rates and mechanisms)
# =============================================================================
make_table2_ari_by_method <- function(results) {

  tbl <- results %>%
    dplyr::group_by(data_label, prox_label, clust_label) %>%
    dplyr::summarise(
      mean_ari = round(mean(ari, na.rm = TRUE), 3),
      sd_ari   = round(sd(ari,   na.rm = TRUE), 3),
      .groups  = "drop"
    ) %>%
    dplyr::mutate(cell = paste0(mean_ari, " (", sd_ari, ")")) %>%
    dplyr::select(data_label, prox_label, clust_label, cell) %>%
    tidyr::pivot_wider(names_from = clust_label, values_from = cell) %>%
    dplyr::rename(
      "Data Type" = data_label,
      "Proximity Measure" = prox_label
    )

  write.csv(tbl, "output/tables/table2_ari_by_clustering_method.csv",
            row.names = FALSE)
  cat("    Saved: output/tables/table2_ari_by_clustering_method.csv\n")
}


# =============================================================================
# FIGURE 1 — Boxplots of ARI by proximity measure and missingness rate
#            Faceted by missingness mechanism (MCAR vs MAR)
#            Shows full distribution, not just means
# =============================================================================
plot_fig1_ari_boxplot_rate <- function(results) {

  # Exclude the complete-data baseline from this figure for cleaner comparison
  plot_data <- dplyr::filter(results, mechanism != "none")

  for (dtype in levels(droplevels(plot_data$data_label))) {
    data_one <- dplyr::filter(plot_data, data_label == dtype)
    file_stem <- DATA_TYPE_FILE_STEMS[[dtype]]

    p <- ggplot2::ggplot(
      data_one,
      ggplot2::aes(x = prox_label, y = ari, fill = miss_pct)
    ) +
      ggplot2::geom_boxplot(
        outlier.size = 0.35, outlier.alpha = 0.25,
        linewidth = 0.45, color = "grey15"
      ) +
      ggplot2::facet_grid(rows = ggplot2::vars(mechanism),
                          cols = ggplot2::vars(separation),
                          labeller = ggplot2::label_both) +
      ggplot2::scale_fill_grey(start = 0.9, end = 0.35, name = "Missingness Rate") +
      ggplot2::scale_y_continuous(breaks = seq(0, 1, 0.25)) +
      ggplot2::coord_cartesian(ylim = c(-0.05, 1.05)) +
      ggplot2::labs(
        x     = "Proximity Measure",
        y     = "Adjusted Rand Index (ARI)"
      ) +
      FIGURE_THEME(base_size = 12) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 30, hjust = 1, size = 10)
      )

    out_file <- paste0("output/figures/fig1_", file_stem, "_ari_boxplot_rate.png")
    ggplot2::ggsave(out_file, p, width = 10, height = 7, dpi = 300)
    cat("    Saved:", out_file, "\n")
  }
}


# =============================================================================
# FIGURE 2 — Faceted dot plot: mean ARI vs. missingness rate.
#            Each proximity measure gets its own panel to avoid crowded
#            black-and-white line encodings.
# =============================================================================
plot_fig2_ari_lineplot_rate <- function(results) {

  # Include the 0% baseline to anchor the left end of each line
  baseline_data <- results %>%
    dplyr::filter(mechanism == "none") %>%
    dplyr::select(-mechanism) %>%
    tidyr::crossing(mechanism = c("MCAR", "MAR"))

  plot_data <- results %>%
    dplyr::filter(mechanism != "none") %>%
    dplyr::bind_rows(baseline_data) %>%
    dplyr::group_by(data_label, prox_label, clust_label, miss_rate, mechanism, separation) %>%
    dplyr::summarise(mean_ari = mean(ari, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(
      mechanism_label = factor(
        paste("Mechanism:", mechanism),
        levels = c("Mechanism: MCAR", "Mechanism: MAR")
      ),
      separation_label = factor(
        paste("Separation:", separation),
        levels = c("Separation: moderate", "Separation: strong")
      )
    )

  for (dtype in levels(droplevels(plot_data$data_label))) {
    data_one <- dplyr::filter(plot_data, data_label == dtype)
    file_stem <- DATA_TYPE_FILE_STEMS[[dtype]]

    p <- ggplot2::ggplot(
      data_one,
      ggplot2::aes(
        x = factor(scales::percent(miss_rate, accuracy = 1),
                   levels = c("0%", "10%", "30%", "50%")),
        y = mean_ari,
        shape = clust_label
      )
    ) +
      ggplot2::geom_point(
        ggplot2::aes(fill = clust_label),
        size = 2.7,
        color = "grey10",
        position = ggplot2::position_dodge(width = 0.45)
      ) +
      ggplot2::facet_grid(
        rows = ggplot2::vars(prox_label),
        cols = ggplot2::vars(mechanism_label, separation_label),
        labeller = ggplot2::label_value,
        switch = "y"
      ) +
      ggplot2::scale_y_continuous(breaks = seq(0, 1, 0.25)) +
      ggplot2::coord_cartesian(ylim = c(-0.05, 1.05)) +
      ggplot2::scale_shape_manual(
        values = c(PAM = 21, UPGMA = 24),
        name   = "Clustering Method"
      ) +
      ggplot2::scale_fill_manual(
        values = c(PAM = "grey15", UPGMA = "white"),
        name   = "Clustering Method"
      ) +
      ggplot2::labs(
        x     = "Missingness Rate",
        y     = "Mean ARI"
      ) +
      FIGURE_THEME(base_size = 11) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(size = 8),
        strip.text.x = ggplot2::element_text(size = 8),
        strip.text.y.left = ggplot2::element_text(angle = 0, size = 8, hjust = 0.5),
        strip.placement = "outside",
        panel.spacing = grid::unit(0.45, "lines")
      )

    out_file <- paste0("output/figures/fig2_", file_stem, "_ari_lineplot_rate.png")
    ggplot2::ggsave(out_file, p, width = 12, height = 8.5, dpi = 300)
    cat("    Saved:", out_file, "\n")
  }
}


# =============================================================================
# FIGURE 3 — Heatmap of mean ARI: rows = proximity measures,
#            columns = missingness rate, panels = mechanism × separation.
#            Cell text shows the mean ARI value.
# =============================================================================
plot_fig3_ari_heatmap <- function(results) {

  heat_data <- results %>%
    dplyr::filter(mechanism != "none") %>%
    dplyr::group_by(data_label, prox_label, miss_pct, mechanism, separation) %>%
    dplyr::summarise(mean_ari = mean(ari, na.rm = TRUE), .groups = "drop")

  for (dtype in levels(droplevels(heat_data$data_label))) {
    data_one <- dplyr::filter(heat_data, data_label == dtype)
    file_stem <- DATA_TYPE_FILE_STEMS[[dtype]]

    p <- ggplot2::ggplot(
      data_one,
      ggplot2::aes(x = miss_pct, y = prox_label, fill = mean_ari)
    ) +
      ggplot2::geom_tile(color = "grey70", linewidth = 0.45) +
      ggplot2::geom_text(ggplot2::aes(label = round(mean_ari, 2)), size = 3.4) +
      ggplot2::facet_grid(rows = ggplot2::vars(mechanism),
                          cols = ggplot2::vars(separation),
                          labeller = ggplot2::label_both) +
      ggplot2::scale_fill_gradient2(
        low      = "white",
        mid      = "grey70",
        high     = "grey20",
        midpoint = 0.5,
        limits   = c(-0.05, 1),
        oob      = scales::squish,
        name     = "Mean ARI"
      ) +
      ggplot2::labs(
        x     = "Missingness Rate",
        y     = "Proximity Measure"
      ) +
      FIGURE_THEME(base_size = 12)

    out_file <- paste0("output/figures/fig3_", file_stem, "_ari_heatmap.png")
    ggplot2::ggsave(out_file, p, width = 9, height = 5.5, dpi = 300)
    cat("    Saved:", out_file, "\n")
  }
}


# =============================================================================
# FIGURE 4 — Silhouette width distributions by proximity measure
#            (secondary metric; faceted by clustering method)
# =============================================================================
plot_fig4_silhouette_boxplot <- function(results) {

  plot_data <- dplyr::filter(results, mechanism != "none")

  for (dtype in levels(droplevels(plot_data$data_label))) {
    data_one <- dplyr::filter(plot_data, data_label == dtype)
    file_stem <- DATA_TYPE_FILE_STEMS[[dtype]]

    p <- ggplot2::ggplot(
      data_one,
      ggplot2::aes(x = prox_label, y = silhouette, fill = miss_pct)
    ) +
      ggplot2::geom_boxplot(
        outlier.size = 0.35, outlier.alpha = 0.25,
        linewidth = 0.45, color = "grey15"
      ) +
      ggplot2::facet_wrap(~ clust_label) +
      ggplot2::scale_fill_grey(start = 0.9, end = 0.35, name = "Missingness Rate") +
      ggplot2::scale_y_continuous(breaks = seq(0, 1, 0.25)) +
      ggplot2::coord_cartesian(ylim = c(-0.1, 1)) +
      ggplot2::labs(
        x     = "Proximity Measure",
        y     = "Mean Silhouette Width"
      ) +
      FIGURE_THEME(base_size = 12) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 30, hjust = 1, size = 10)
      )

    out_file <- paste0("output/figures/fig4_", file_stem, "_silhouette_boxplot.png")
    ggplot2::ggsave(out_file, p, width = 9, height = 5.5, dpi = 300)
    cat("    Saved:", out_file, "\n")
  }
}
