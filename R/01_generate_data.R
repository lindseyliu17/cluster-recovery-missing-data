# =============================================================================
# 01_generate_data.R
#
# PURPOSE: Generate synthetic data with a known three-cluster structure.
#          Because the data are simulated, we know the TRUE cluster label for
#          every observation — this lets us measure later how well the
#          clustering algorithms recover that truth.
#
# KEY IDEA: The continuous condition uses a multivariate-normal mixture.  The
#           mixed condition starts from the same latent structure, then creates
#           continuous, ordinal, binary, and nominal observed variables.
#
#           In both conditions, an auxiliary fully observed MAR driver is
#           generated as a noisy linear function of the first cluster-informative
#           latent dimension. It is NOT used as a clustering variable. Because
#           the driver correlates with cluster membership, MAR missingness is
#           systematically heavier in some clusters than others — making MAR
#           meaningfully harder than MCAR rather than equivalent to it.
# =============================================================================


#' Generate data from a three-cluster simulation model
#'
#' @param n_per_cluster  Number of observations per cluster (total n = 3 × this)
#' @param n_vars         Number of clustering variables. Continuous data use this
#'                       value directly; mixed data currently use 5 variables.
#' @param separation     "moderate" (delta = 2 SD) or "strong" (delta = 4 SD)
#' @param data_type      "continuous" or "mixed"
#'
#' @return A list with two elements:
#'   $data       — clustering variables, rows = observations
#'   $labels     — integer vector of true cluster memberships (1, 2, or 3)
#'   $mar_driver — fully observed auxiliary variable for MAR missingness
generate_cluster_data <- function(n_per_cluster, n_vars, separation,
                                  data_type = "continuous") {

  # --- Input validation -------------------------------------------------------
  if (!separation %in% c("moderate", "strong")) {
    stop("separation must be 'moderate' or 'strong'")
  }
  if (!data_type %in% c("continuous", "mixed")) {
    stop("data_type must be 'continuous' or 'mixed'")
  }

  # --- Cluster separation distance --------------------------------------------
  # delta = how many SD units separate adjacent cluster centers.
  # Moderate (delta=2): clusters overlap noticeably — harder recovery task.
  # Strong  (delta=4): clusters are well-separated — easier recovery task.
  delta <- switch(separation,
    moderate = 2.0,
    strong   = 4.0
  )

  # --- Cluster means ----------------------------------------------------------
  # Centers are arranged in an equilateral triangle on dimensions 1 and 2.
  # Dimensions 3:p have all cluster means at zero (no separation there).
  #
  #   Cluster 1 at (0, 0)
  #   Cluster 2 at (delta, 0)
  #   Cluster 3 at (delta/2,  delta * sqrt(3)/2)
  #
  # This gives each pair of clusters the same Euclidean distance (= delta).
  cluster_means <- matrix(0, nrow = 3, ncol = n_vars)
  cluster_means[1, 1:2] <- c(0,         0)
  cluster_means[2, 1:2] <- c(delta,     0)
  cluster_means[3, 1:2] <- c(delta / 2, delta * sqrt(3) / 2)

  # --- Covariance matrix ------------------------------------------------------
  # Identity matrix: each variable has variance 1, no cross-variable correlations.
  sigma <- diag(n_vars)

  # --- Simulate latent continuous data ----------------------------------------
  # Draw n_per_cluster observations from each cluster's normal distribution.
  cluster1 <- MASS::mvrnorm(n = n_per_cluster, mu = cluster_means[1, ], Sigma = sigma)
  cluster2 <- MASS::mvrnorm(n = n_per_cluster, mu = cluster_means[2, ], Sigma = sigma)
  cluster3 <- MASS::mvrnorm(n = n_per_cluster, mu = cluster_means[3, ], Sigma = sigma)

  # Stack rows and create the matching true-label vector
  latent <- rbind(cluster1, cluster2, cluster3)
  labels <- rep(1:3, each = n_per_cluster)

  # --- Convert latent variables into the requested observed data type ----------
  if (data_type == "continuous") {
    data <- latent
    colnames(data) <- paste0("x", seq_len(ncol(data)))
  } else {
    data <- make_mixed_cluster_data(latent, labels)
  }

  # Fully observed auxiliary variable used only to generate MAR missingness.
  # Defined as the first cluster-informative latent dimension plus noise, so it
  # is correlated with cluster membership but is not an analysis variable itself.
  # Using a cluster-correlated driver makes MAR missingness heavier in some
  # clusters than others — genuinely distinguishable from MCAR.
  mar_driver <- latent[, 1] + rnorm(nrow(latent))

  # Shuffle rows so clusters are not in contiguous blocks.
  # Some algorithms are sensitive to data ordering; shuffling avoids that.
  shuffle_idx <- sample(nrow(latent))

  return(list(
    data       = data[shuffle_idx, ],
    labels     = labels[shuffle_idx],
    mar_driver = mar_driver[shuffle_idx]
  ))
}


# =============================================================================
# Mixed-type observed data
#
# Variables:
#   x_cont1, x_cont2  — continuous cluster-informative variables
#   x_ord             — ordered three-level variable derived from latent signal
#   x_bin             — binary factor with cluster-varying probabilities
#   x_nom             — nominal factor with cluster-varying probabilities
# =============================================================================
make_mixed_cluster_data <- function(latent, labels) {

  n <- nrow(latent)

  x_cont1 <- latent[, 1]
  x_cont2 <- latent[, 2]

  # Ordered categorical variable: cluster signal plus noise, then tertiles.
  ord_score <- 0.65 * latent[, 1] + 0.35 * latent[, 2] + rnorm(n, sd = 0.75)
  ord_breaks <- stats::quantile(ord_score, probs = c(0, 1/3, 2/3, 1))
  ord_breaks <- unique(ord_breaks)
  if (length(ord_breaks) < 4) {
    ord_breaks <- seq(min(ord_score), max(ord_score), length.out = 4)
  }
  x_ord <- cut(
    ord_score,
    breaks = ord_breaks,
    include.lowest = TRUE,
    labels = c("low", "medium", "high"),
    ordered_result = TRUE
  )

  # Binary categorical variable with different probabilities by true cluster.
  bin_probs <- c(0.20, 0.50, 0.80)[labels]
  x_bin <- factor(
    ifelse(runif(n) < bin_probs, "yes", "no"),
    levels = c("no", "yes")
  )

  # Nominal categorical variable with cluster-specific category profiles.
  nom_prob_mat <- rbind(
    c(0.70, 0.20, 0.10),
    c(0.20, 0.60, 0.20),
    c(0.10, 0.25, 0.65)
  )
  x_nom <- character(n)
  for (i in seq_len(n)) {
    x_nom[i] <- sample(c("A", "B", "C"), size = 1, prob = nom_prob_mat[labels[i], ])
  }
  x_nom <- factor(x_nom, levels = c("A", "B", "C"))

  data.frame(
    x_cont1 = x_cont1,
    x_cont2 = x_cont2,
    x_ord   = x_ord,
    x_bin   = x_bin,
    x_nom   = x_nom
  )
}
