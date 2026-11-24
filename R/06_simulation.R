# =============================================================================
# 06_simulation.R
#
# PURPOSE: Orchestrate a single simulation replicate and the full simulation
#          loop across all experimental conditions.
#
# FLOW OF ONE REPLICATE:
#   1. Generate complete cluster data                     (01_generate_data.R)
#   2. Introduce missing values                            (02_missingness.R)
#   3. Compute 5 proximity matrices                        (03_proximity.R)
#   4. Apply 2 clustering methods to each matrix           (04_clustering.R)
#   5. Evaluate ARI and silhouette for all 10 combinations (05_evaluation.R)
#   → Return a 10-row data frame (one row per proximity × clustering combo)
#
# The full simulation loops over all conditions × replicates and stacks
# results into one large data frame that is saved to disk.
# =============================================================================


#' Run one simulation replicate for a given set of conditions
#'
#' @param separation     Cluster separation: "moderate" or "strong"
#' @param mechanism      Missingness mechanism: "MCAR", "MAR", or "none"
#' @param miss_rate      Proportion of missing values (0, 0.10, 0.30, 0.50)
#' @param n_per_cluster  Observations per cluster
#' @param n_vars         Number of variables
#' @param n_clusters     Number of clusters (= 3, the true value)
#' @param data_type      "continuous" or "mixed"
#'
#' @return Data frame with 10 rows (5 proximity measures × 2 clustering methods)
#'         and columns: proximity, clustering, ari, silhouette
run_one_replicate <- function(separation, mechanism, miss_rate,
                               n_per_cluster, n_vars, n_clusters,
                               data_type = "continuous") {

  # Step 1 — Generate complete data with known cluster labels
  sim_data    <- generate_cluster_data(n_per_cluster, n_vars, separation, data_type)
  X_complete  <- sim_data$data
  true_labels <- sim_data$labels

  # Step 2 — Introduce missing values (returns data unchanged when miss_rate = 0)
  X_missing <- introduce_missingness(
    X_complete,
    mechanism,
    miss_rate,
    mar_driver = sim_data$mar_driver
  )

  # Step 3 — Compute all five proximity matrices from the (possibly incomplete) data
  prox_list <- compute_all_proximity(X_missing)

  # Step 4 & 5 — Cluster with each method and evaluate each combination
  proximity_names <- names(prox_list)
  cluster_methods <- c("pam", "average")

  # Pre-allocate a list to hold results (faster than growing a data frame in a loop)
  n_combos <- length(proximity_names) * length(cluster_methods)
  result_list <- vector("list", n_combos)
  row_idx <- 1

  for (prox_name in proximity_names) {
    dist_mat <- prox_list[[prox_name]]

    # Number/proportion of observation pairs with no jointly observed variable
    # (defined for pairwise-complete methods and Gower; NA for imputed methods).
    zero_pairs <- attr(dist_mat, "n_zero_overlap_pairs")
    if (is.null(zero_pairs)) zero_pairs <- NA_integer_
    zero_pair_prop <- attr(dist_mat, "prop_zero_overlap_pairs")
    if (is.null(zero_pair_prop)) zero_pair_prop <- NA_real_

    for (clust_method in cluster_methods) {

      # Wrap in tryCatch so a single failed replicate does not crash the whole simulation.
      # If something goes wrong (e.g., PAM fails to converge), we record NA for metrics.
      metrics <- tryCatch({
        est_labels <- cluster_from_distance(dist_mat, clust_method, k = n_clusters)
        list(
          ari        = compute_ari(true_labels, est_labels),
          silhouette = compute_silhouette(est_labels, dist_mat)
        )
      }, error = function(e) {
        list(ari = NA_real_, silhouette = NA_real_)
      })

      result_list[[row_idx]] <- data.frame(
        proximity  = prox_name,
        clustering = clust_method,
        ari        = metrics$ari,
        silhouette = metrics$silhouette,
        zero_pairs = zero_pairs,
        zero_pair_prop = zero_pair_prop,
        stringsAsFactors = FALSE
      )
      row_idx <- row_idx + 1
    }
  }

  do.call(rbind, result_list)
}


# =============================================================================

#' Run the complete simulation study across all conditions and replicates
#'
#' @param conditions     Data frame of experimental conditions (from expand.grid)
#'                       Must have columns: separation, mechanism, miss_rate
#' @param n_sim          Number of simulation replicates per condition
#' @param n_per_cluster  Observations per cluster
#' @param n_vars         Number of variables
#' @param n_clusters     Number of clusters
#'
#' @return Data frame with all results; columns include condition labels and
#'         replicate index alongside ari and silhouette
run_full_simulation <- function(conditions, n_sim,
                                 n_per_cluster, n_vars, n_clusters) {

  n_conditions <- nrow(conditions)
  total_reps   <- n_conditions * n_sim

  # Pre-allocate outer list: one entry per condition × replicate combination
  all_results <- vector("list", total_reps)
  global_idx  <- 1

  for (cond_i in seq_len(n_conditions)) {

    sep   <- conditions$separation[cond_i]
    mech  <- conditions$mechanism[cond_i]
    mrate <- conditions$miss_rate[cond_i]
    dtype <- if ("data_type" %in% names(conditions)) {
      conditions$data_type[cond_i]
    } else {
      "continuous"
    }

    cat(sprintf(
      "[Condition %2d / %2d]  data = %-10s | separation = %-8s | mechanism = %-4s | miss_rate = %.0f%%\n",
      cond_i, n_conditions, dtype, sep, mech, mrate * 100
    ))

    for (rep_i in seq_len(n_sim)) {

      # Lightweight progress tick every 50 replicates within a condition
      if (rep_i %% 50 == 0) {
        cat(sprintf("  ... replicate %d / %d\n", rep_i, n_sim))
      }

      rep_result <- run_one_replicate(
        separation    = sep,
        mechanism     = mech,
        miss_rate     = mrate,
        n_per_cluster = n_per_cluster,
        n_vars        = n_vars,
        n_clusters    = n_clusters,
        data_type     = dtype
      )

      # Tag each row with the condition labels and replicate number
      rep_result$separation <- sep
      rep_result$mechanism  <- mech
      rep_result$miss_rate  <- mrate
      rep_result$data_type  <- dtype
      rep_result$replicate  <- rep_i

      all_results[[global_idx]] <- rep_result
      global_idx <- global_idx + 1
    }
  }

  cat("\nSimulation complete.\n")
  do.call(rbind, all_results)
}
