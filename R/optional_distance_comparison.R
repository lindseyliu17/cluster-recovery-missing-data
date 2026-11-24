# =============================================================================
# optional_distance_comparison.R
#
# PURPOSE: Compare distance matrices computed from COMPLETE data vs. the same
#          matrices recomputed after introducing missing values.  This shows
#          HOW MUCH each proximity measure distorts the underlying similarity
#          structure when data are incomplete.
#
# NOT PART OF THE MAIN WORKFLOW — run independently:
#
#   source("R/01_generate_data.R")
#   source("R/02_missingness.R")
#   source("R/03_proximity.R")
#   source("R/optional_distance_comparison.R")
#   run_distance_comparison()          # uses default settings
#
# OUTPUT (saved to output/optional/):
#   distance_comparison_metrics.csv   — numeric summary (correlation, RMSE, etc.)
#   fig_dist_scatter_<tag>.png        — scatter plots: complete vs incomplete
#   fig_dist_metrics_summary.png      — faceted line plot across conditions
# =============================================================================


# =============================================================================
# METRIC FUNCTIONS
# =============================================================================

#' Extract the upper triangle of a distance matrix as a plain numeric vector
#'
#' Both matrices must have the same dimensions.  Diagonal is excluded because
#' all self-distances are zero by definition and carry no information.
upper_tri_vec <- function(d_mat) {
  d_mat[upper.tri(d_mat)]
}


#' Compute pairwise comparison metrics between two distance matrices
#'
#' @param d_complete   Full-data distance matrix (as dist or matrix)
#' @param d_incomplete Incomplete-data distance matrix (same dimensions)
#'
#' @return Named numeric vector with:
#'   pearson_r   — Pearson correlation of the upper-triangle entries
#'   spearman_r  — Spearman rank correlation (rank-order preservation)
#'   rmse        — Root mean squared error (absolute distortion)
#'   rel_rmse    — RMSE relative to the mean complete-data distance (scale-free)
#'   mae         — Mean absolute error
compare_distance_matrices <- function(d_complete, d_incomplete) {

  # Convert to plain matrices so indexing is consistent
  mat_c <- as.matrix(d_complete)
  mat_i <- as.matrix(d_incomplete)

  if (!identical(dim(mat_c), dim(mat_i))) {
    stop("Distance matrices must have the same dimensions")
  }

  v_c <- upper_tri_vec(mat_c)   # complete-data distances, upper triangle
  v_i <- upper_tri_vec(mat_i)   # incomplete-data distances, upper triangle

  # Keep only pairs where both values are finite (guards against any residual NAs)
  keep <- is.finite(v_c) & is.finite(v_i)
  v_c  <- v_c[keep]
  v_i  <- v_i[keep]

  mean_c  <- mean(v_c)
  errors  <- v_i - v_c

  c(
    pearson_r  = cor(v_c, v_i, method = "pearson"),
    spearman_r = cor(v_c, v_i, method = "spearman"),
    rmse       = sqrt(mean(errors^2)),
    rel_rmse   = sqrt(mean(errors^2)) / mean_c,
    mae        = mean(abs(errors))
  )
}


# =============================================================================
# SCATTER PLOT: complete vs incomplete distances for one matrix pair
# =============================================================================

#' Scatter plot comparing complete-data vs incomplete-data distances
#'
#' Each point is one (i, j) pair from the upper triangle of the distance matrix.
#' The grey diagonal is the "no distortion" reference line.
#'
#' @param d_complete     Complete-data distance matrix
#' @param d_incomplete   Incomplete-data distance matrix
#' @param title          Plot title string
#' @param n_sample       Max number of pairs to plot (random sample if exceeded)
#'                       Keeps the scatter readable when n is large
#'
#' @return A ggplot2 object
plot_distance_scatter <- function(d_complete, d_incomplete,
                                  title = "Complete vs Incomplete Distances",
                                  n_sample = 3000) {

  mat_c <- as.matrix(d_complete)
  mat_i <- as.matrix(d_incomplete)

  v_c <- upper_tri_vec(mat_c)
  v_i <- upper_tri_vec(mat_i)

  keep <- is.finite(v_c) & is.finite(v_i)
  df   <- data.frame(complete = v_c[keep], incomplete = v_i[keep])

  # Randomly down-sample if there are too many points to display clearly
  if (nrow(df) > n_sample) {
    df <- df[sample(nrow(df), n_sample), ]
  }

  pearson_label <- sprintf("Pearson r = %.3f", cor(df$complete, df$incomplete))

  ggplot2::ggplot(df, ggplot2::aes(x = complete, y = incomplete)) +
    ggplot2::geom_point(alpha = 0.15, size = 0.6, color = "#2166ac") +
    ggplot2::geom_abline(slope = 1, intercept = 0,
                         color = "grey40", linewidth = 0.6, linetype = "dashed") +
    ggplot2::annotate("text",
                      x = Inf, y = -Inf, hjust = 1.05, vjust = -0.5,
                      label = pearson_label, size = 3.2, color = "grey30") +
    ggplot2::labs(
      title    = title,
      subtitle = "Dashed line = perfect agreement (no distortion)",
      x        = "Complete-data distance",
      y        = "Incomplete-data distance"
    ) +
    ggplot2::theme_bw(base_size = 10)
}


# =============================================================================
# SUMMARY LINE PLOT: metrics across proximity measures and missingness rates
# =============================================================================

#' Line plot of Pearson r and relative RMSE across missingness rates
#'
#' @param metrics_df  Data frame returned by run_distance_comparison()
plot_distance_metrics_summary <- function(metrics_df) {

  # Pivot to long format for ggplot faceting
  long_df <- metrics_df %>%
    dplyr::select(proximity, mechanism, miss_rate, pearson_r, rel_rmse) %>%
    tidyr::pivot_longer(
      cols      = c("pearson_r", "rel_rmse"),
      names_to  = "metric",
      values_to = "value"
    ) %>%
    dplyr::mutate(
      metric = dplyr::recode(metric,
        pearson_r = "Pearson r (higher = better)",
        rel_rmse  = "Relative RMSE (lower = better)"
      ),
      prox_label = dplyr::recode(proximity,
        euclidean_imputed  = "Eucl. (Imputed)",
        euclidean_pairwise = "Eucl. (Pairwise)",
        manhattan_imputed  = "Manh. (Imputed)",
        manhattan_pairwise = "Manh. (Pairwise)",
        gower              = "Gower"
      ),
      prox_label = factor(prox_label, levels = c(
        "Eucl. (Imputed)", "Eucl. (Pairwise)",
        "Manh. (Imputed)", "Manh. (Pairwise)",
        "Gower"
      ))
    )

  ggplot2::ggplot(
    long_df,
    ggplot2::aes(x = miss_rate, y = value,
                 color = prox_label, group = prox_label, linetype = prox_label)
  ) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 2) +
    ggplot2::facet_grid(
      rows     = ggplot2::vars(metric),
      cols     = ggplot2::vars(mechanism),
      scales   = "free_y"
    ) +
    ggplot2::scale_x_continuous(labels = scales::percent,
                                breaks = c(0.10, 0.30, 0.50)) +
    ggplot2::scale_color_brewer(palette = "Dark2", name = "Proximity Measure") +
    ggplot2::scale_linetype_manual(
      values = c("solid", "dashed", "dotdash", "dotted", "longdash"),
      name   = "Proximity Measure"
    ) +
    ggplot2::labs(
      title = "How Much Does Missingness Distort the Distance Matrix?",
      x     = "Missingness Rate",
      y     = "Metric Value"
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(legend.position = "bottom")
}


# =============================================================================
# ORCHESTRATOR
# =============================================================================

#' Run the full distance-matrix comparison across conditions
#'
#' For each combination of (data_type × separation × mechanism × miss_rate),
#' generates one dataset, computes the COMPLETE-data distance matrices, then
#' introduces missingness and recomputes them.  Numeric metrics are saved to
#' CSV; scatter plots (one per proximity × condition) and a summary figure are
#' saved to output/optional/.
#'
#' @param n_per_cluster   Observations per cluster (default 100; total n = 300)
#' @param n_vars          Number of clustering variables (default 5)
#' @param seed            Random seed for reproducibility
#' @param data_types      Character vector: "continuous", "mixed", or both
#' @param separations     Character vector: "moderate", "strong", or both
#' @param mechanisms      Character vector: "MCAR", "MAR", or both
#' @param miss_rates      Numeric vector of missingness proportions
#' @param save_scatter    Logical — save individual scatter PNGs (can be many)
run_distance_comparison <- function(
    n_per_cluster = 100,
    n_vars        = 5,
    seed          = 2024,
    data_types    = c("continuous", "mixed"),
    separations   = c("moderate", "strong"),
    mechanisms    = c("MCAR", "MAR"),
    miss_rates    = c(0.10, 0.30, 0.50),
    save_scatter  = TRUE
) {

  dir.create("output/optional", recursive = TRUE, showWarnings = FALSE)

  set.seed(seed)

  # Build the condition grid
  conditions <- expand.grid(
    data_type  = data_types,
    separation = separations,
    mechanism  = mechanisms,
    miss_rate  = miss_rates,
    stringsAsFactors = FALSE
  )

  # Storage for the summary metrics (one row per condition × proximity measure)
  all_metrics <- vector("list", nrow(conditions) * 5)
  row_idx     <- 1

  cat(sprintf("Running distance comparison across %d conditions...\n", nrow(conditions)))

  for (cond_i in seq_len(nrow(conditions))) {

    dtype  <- conditions$data_type[cond_i]
    sep    <- conditions$separation[cond_i]
    mech   <- conditions$mechanism[cond_i]
    mrate  <- conditions$miss_rate[cond_i]

    cat(sprintf(
      "  [%2d/%2d]  %-10s | %-8s | %-4s | %.0f%%\n",
      cond_i, nrow(conditions), dtype, sep, mech, mrate * 100
    ))

    # Step 1 — Generate one complete dataset
    sim_data    <- generate_cluster_data(n_per_cluster, n_vars, sep, dtype)
    X_complete  <- sim_data$data

    # Step 2 — Distance matrices from COMPLETE data (ground truth)
    prox_complete <- compute_all_proximity(X_complete)

    # Step 3 — Introduce missingness
    X_missing <- introduce_missingness(
      X_complete, mech, mrate,
      mar_driver = sim_data$mar_driver
    )

    # Step 4 — Distance matrices from INCOMPLETE data
    prox_missing <- compute_all_proximity(X_missing)

    # Step 5 — Compare each proximity measure
    for (prox_name in names(prox_complete)) {

      metrics_vec <- compare_distance_matrices(
        prox_complete[[prox_name]],
        prox_missing[[prox_name]]
      )

      all_metrics[[row_idx]] <- data.frame(
        data_type  = dtype,
        separation = sep,
        mechanism  = mech,
        miss_rate  = mrate,
        proximity  = prox_name,
        as.data.frame(t(metrics_vec)),
        stringsAsFactors = FALSE
      )
      row_idx <- row_idx + 1

      # Step 6 — (Optional) scatter plot for this specific condition
      if (save_scatter) {

        scatter_title <- sprintf(
          "%s | %s | %s | %.0f%% missing\n%s distance: complete vs incomplete",
          dtype, sep, mech, mrate * 100, prox_name
        )

        p_scatter <- plot_distance_scatter(
          prox_complete[[prox_name]],
          prox_missing[[prox_name]],
          title = scatter_title
        )

        scatter_file <- sprintf(
          "output/optional/dist_scatter_%s_%s_%s_%02.0fpct_%s.png",
          dtype, sep, tolower(mech), mrate * 100, prox_name
        )
        ggplot2::ggsave(scatter_file, p_scatter, width = 5, height = 4, dpi = 120)
      }
    }
  }

  # Combine all metrics into one data frame
  metrics_df <- do.call(rbind, all_metrics)

  # Save numeric summary
  write.csv(metrics_df, "output/optional/distance_comparison_metrics.csv",
            row.names = FALSE)
  cat("\n  Saved: output/optional/distance_comparison_metrics.csv\n")

  # Save summary figure
  p_summary <- plot_distance_metrics_summary(metrics_df)
  ggplot2::ggsave(
    "output/optional/fig_dist_metrics_summary.png",
    p_summary, width = 10, height = 6, dpi = 150
  )
  cat("  Saved: output/optional/fig_dist_metrics_summary.png\n")

  if (save_scatter) {
    cat(sprintf(
      "  Scatter plots saved to output/optional/ (%d files)\n",
      nrow(conditions) * length(names(prox_complete))
    ))
  }

  cat("\nDone.\n")
  invisible(metrics_df)
}
