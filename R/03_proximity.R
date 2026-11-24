# =============================================================================
# 03_proximity.R
#
# PURPOSE: Compute pairwise distance / dissimilarity matrices from data that
#          may contain missing values.  Five measures are implemented:
#
#   1. euclidean_imputed   — replace NAs with column means, then Euclidean
#   2. euclidean_pairwise  — Euclidean using only variables observed in BOTH
#                            members of each pair (pairwise-complete approach)
#   3. manhattan_imputed   — mean imputation, then Manhattan (L1) distance
#   4. manhattan_pairwise  — Manhattan using pairwise-complete variables
#   5. gower               — Gower dissimilarity, which handles NAs natively
#                            by skipping missing components per pair
#
# MIXED-TYPE NOTE:
#   For mixed data, Euclidean and Manhattan distances are computed after a
#   simple numeric encoding: continuous variables stay numeric, ordered factors
#   become integer scores, and unordered factors are expanded into dummy
#   indicators. These measures are included as comparison methods. Gower is the
#   proximity measure designed to handle mixed continuous/categorical data
#   directly.
#
# SCALING NOTE (pairwise methods):
#   When fewer than p variables are available for a pair, a raw sum of
#   differences would produce a smaller distance than if all p variables
#   contributed — even if the clusters are equally separated.  To keep
#   distances on the same scale regardless of how many variables are missing,
#   we multiply by p/k (Euclidean: inside the sqrt via p/k factor on the sum
#   of squares; Manhattan: direct multiplication).  This is sometimes called
#   the "normalized pairwise" approach.
#
# NA HANDLING (pairwise methods):
#   If two observations share NO jointly-observed variable (k = 0), the
#   distance is undefined.  We replace these undefined entries with the
#   maximum finite distance in the matrix — a conservative assumption
#   (treat unknowns as far apart).
#
# All functions return an object of class "dist" (lower-triangular matrix),
# which is the input format expected by cluster::pam() and stats::hclust().
# =============================================================================


#' Compute all five proximity matrices for a dataset with missing values
#'
#' @param data  Numeric matrix (may contain NAs)
#' @return Named list of five "dist" objects
compute_all_proximity <- function(data) {
  list(
    euclidean_imputed  = proximity_euclidean_imputed(data),
    euclidean_pairwise = proximity_euclidean_pairwise(data),
    manhattan_imputed  = proximity_manhattan_imputed(data),
    manhattan_pairwise = proximity_manhattan_pairwise(data),
    gower              = proximity_gower(data)
  )
}


# =============================================================================
# Internal helper: mean imputation
# Replace every NA with the mean of non-missing values in that column.
# =============================================================================
impute_column_means <- function(data) {
  data_imp <- data
  for (j in seq_len(ncol(data))) {
    col_mean <- mean(data[, j], na.rm = TRUE)
    if (is.nan(col_mean)) col_mean <- 0
    data_imp[is.na(data_imp[, j]), j] <- col_mean
  }
  return(data_imp)
}


# =============================================================================
# Internal helper: standardize encoded columns to mean 0, SD 1
#
# Computed from observed (non-NA) values per column so NAs are not affected.
# Called after encoding and before any distance calculation so that continuous
# variables (spanning several units) and binary dummy indicators (0/1) enter
# distances on a comparable scale.  Constant or near-constant columns are
# skipped to avoid division by zero.
# =============================================================================
standardize_encoded_columns <- function(data) {
  for (j in seq_len(ncol(data))) {
    obs_vals <- data[!is.na(data[, j]), j]
    if (length(obs_vals) < 2) next
    s <- sd(obs_vals)
    if (s < 1e-10) next
    data[, j] <- (data[, j] - mean(obs_vals)) / s
  }
  data
}


# =============================================================================
# Internal helper: replace undefined pair distances
#
# When a pair has no jointly observed variables, its distance is undefined.
# Replace undefined off-diagonal entries with the largest finite off-diagonal
# distance. If no finite off-diagonal distance exists, use 1 as a neutral
# positive fallback rather than using the diagonal zeros.
# =============================================================================
replace_undefined_distances <- function(d_mat) {
  off_diag <- upper.tri(d_mat) | lower.tri(d_mat)
  finite_off_diag <- d_mat[off_diag & is.finite(d_mat)]
  max_d <- if (length(finite_off_diag) > 0) max(finite_off_diag) else 1

  d_mat[is.na(d_mat)] <- max_d
  diag(d_mat) <- 0
  d_mat
}


# =============================================================================
# Internal helper: numeric encoding for Euclidean / Manhattan distances
# =============================================================================
encode_for_numeric_distance <- function(data) {
  data_df <- as.data.frame(data)

  encoded_cols <- list()
  for (j in seq_along(data_df)) {
    x <- data_df[[j]]
    var_name <- names(data_df)[j]
    if (is.null(var_name) || var_name == "") var_name <- paste0("x", j)

    if (is.ordered(x)) {
      encoded_cols[[var_name]] <- as.numeric(x)
    } else if (is.factor(x) || is.character(x) || is.logical(x)) {
      f <- factor(x)
      for (lev in levels(f)) {
        dummy_name <- paste(var_name, make.names(lev), sep = "_")
        encoded_cols[[dummy_name]] <- ifelse(
          is.na(f),
          NA_real_,
          as.numeric(f == lev)
        )
      }
    } else if (is.numeric(x) || is.integer(x)) {
      encoded_cols[[var_name]] <- as.numeric(x)
    } else {
      f <- factor(x)
      for (lev in levels(f)) {
        dummy_name <- paste(var_name, make.names(lev), sep = "_")
        encoded_cols[[dummy_name]] <- ifelse(
          is.na(f),
          NA_real_,
          as.numeric(f == lev)
        )
      }
    }
  }

  as.matrix(as.data.frame(encoded_cols, check.names = FALSE))
}


# =============================================================================
# Internal helper: count upper-triangle pairs with no jointly observed column
# =============================================================================
count_zero_overlap_pairs <- function(data) {
  obs_mat <- as.matrix(!is.na(as.data.frame(data))) * 1L
  k_mat   <- obs_mat %*% t(obs_mat)
  sum(k_mat[upper.tri(k_mat)] == 0L)
}


# =============================================================================
# 1. Euclidean distance after mean imputation
# =============================================================================
proximity_euclidean_imputed <- function(data) {
  data     <- encode_for_numeric_distance(data)
  data     <- standardize_encoded_columns(data)
  data_imp <- impute_column_means(data)
  dist(data_imp, method = "euclidean")
}


# =============================================================================
# 2. Euclidean distance — pairwise-complete approach
#
# For each pair (i, j):
#   k    = number of variables observed in BOTH i and j
#   d_ij = sqrt( (p/k) * sum_{jointly observed} (x_il - x_jl)^2 )
#
# Vectorized implementation: accumulate squared differences variable-by-
# variable, then apply the p/k scaling at the end.
# =============================================================================
proximity_euclidean_pairwise <- function(data) {
  data <- encode_for_numeric_distance(data)
  data <- standardize_encoded_columns(data)
  n <- nrow(data)
  p <- ncol(data)

  sum_sq <- matrix(0, nrow = n, ncol = n)

  for (j in seq_len(p)) {
    obs_j <- which(!is.na(data[, j]))
    if (length(obs_j) < 2) next

    x_j       <- data[obs_j, j]
    diff_sq_j <- outer(x_j, x_j, "-")^2
    sum_sq[obs_j, obs_j] <- sum_sq[obs_j, obs_j] + diff_sq_j
  }

  obs_mat <- (!is.na(data)) * 1L
  k_mat   <- obs_mat %*% t(obs_mat)

  # Count upper-triangle pairs with no jointly observed variable.
  # These distances are undefined; they are later replaced with the max finite
  # distance as a conservative fallback.
  n_zero_pairs <- sum(k_mat[upper.tri(k_mat)] == 0L)

  scale_mat          <- matrix(NA_real_, n, n)
  nonzero            <- k_mat > 0
  scale_mat[nonzero] <- p / k_mat[nonzero]

  d_mat <- sqrt(sum_sq * scale_mat)

  d_mat <- replace_undefined_distances(d_mat)

  d_out <- as.dist(d_mat)
  attr(d_out, "n_zero_overlap_pairs") <- n_zero_pairs
  attr(d_out, "prop_zero_overlap_pairs") <- n_zero_pairs / choose(n, 2)
  d_out
}


# =============================================================================
# 3. Manhattan distance after mean imputation
# =============================================================================
proximity_manhattan_imputed <- function(data) {
  data     <- encode_for_numeric_distance(data)
  data     <- standardize_encoded_columns(data)
  data_imp <- impute_column_means(data)
  dist(data_imp, method = "manhattan")
}


# =============================================================================
# 4. Manhattan distance — pairwise-complete approach
#
# For each pair (i, j):
#   d_ij = (p/k) * sum_{jointly observed} |x_il - x_jl|
#
# Same vectorized structure as the Euclidean pairwise method above.
# =============================================================================
proximity_manhattan_pairwise <- function(data) {
  data <- encode_for_numeric_distance(data)
  data <- standardize_encoded_columns(data)
  n <- nrow(data)
  p <- ncol(data)

  sum_abs <- matrix(0, nrow = n, ncol = n)

  for (j in seq_len(p)) {
    obs_j <- which(!is.na(data[, j]))
    if (length(obs_j) < 2) next

    x_j        <- data[obs_j, j]
    abs_diff_j <- abs(outer(x_j, x_j, "-"))
    sum_abs[obs_j, obs_j] <- sum_abs[obs_j, obs_j] + abs_diff_j
  }

  obs_mat <- (!is.na(data)) * 1L
  k_mat   <- obs_mat %*% t(obs_mat)

  n_zero_pairs <- sum(k_mat[upper.tri(k_mat)] == 0L)

  scale_mat          <- matrix(NA_real_, n, n)
  nonzero            <- k_mat > 0
  scale_mat[nonzero] <- p / k_mat[nonzero]

  d_mat <- sum_abs * scale_mat

  d_mat <- replace_undefined_distances(d_mat)

  d_out <- as.dist(d_mat)
  attr(d_out, "n_zero_overlap_pairs") <- n_zero_pairs
  attr(d_out, "prop_zero_overlap_pairs") <- n_zero_pairs / choose(n, 2)
  d_out
}


# =============================================================================
# 5. Gower dissimilarity
#
# cluster::daisy() with metric = "gower" computes Gower's general similarity
# coefficient.  For each pair (i, j) and each variable, it contributes a
# partial similarity only when BOTH observations are non-missing.  The final
# dissimilarity is 1 − (weighted average of partial similarities).
#
# No imputation is needed — Gower handles the missing-data weighting internally
# when a pair shares at least one observed variable. If a pair shares no observed
# variables, daisy() returns NA; we replace those entries with the maximum finite
# dissimilarity so clustering algorithms can still use the matrix.
# =============================================================================
proximity_gower <- function(data) {
  n <- nrow(data)
  n_zero_pairs <- count_zero_overlap_pairs(data)

  # daisy() expects a data frame; the result is class "dissimilarity"
  # which as.dist() converts to standard "dist".
  diss <- cluster::daisy(as.data.frame(data), metric = "gower")
  d_mat <- as.matrix(as.dist(diss))

  if (anyNA(d_mat)) d_mat <- replace_undefined_distances(d_mat)

  d_out <- as.dist(d_mat)
  attr(d_out, "n_zero_overlap_pairs") <- n_zero_pairs
  attr(d_out, "prop_zero_overlap_pairs") <- n_zero_pairs / choose(n, 2)
  d_out
}


# =============================================================================
# Scaling verification
#
# Call once after sourcing this file to confirm that pairwise Euclidean and
# Manhattan distances equal their standard counterparts on complete data.
# On complete data k_ij = p for every pair, so p/k = 1 and the pairwise
# formula reduces to the ordinary distance.  After standardisation both
# paths operate on identically scaled columns, so results should match to
# floating-point precision.
#
# Usage:
#   source("R/03_proximity.R")
#   verify_pairwise_scaling()
# =============================================================================
verify_pairwise_scaling <- function(n = 50, p = 5, seed = 1) {

  set.seed(seed)
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(X) <- paste0("x", seq_len(p))

  # Both functions standardise internally, so compare against dist() on the
  # same standardised matrix.
  X_std <- standardize_encoded_columns(X)

  check <- function(label, d_pairwise, method) {
    d_ref  <- as.matrix(dist(X_std, method = method))
    d_pair <- as.matrix(d_pairwise)
    max_diff <- max(abs(d_pair - d_ref))
    if (max_diff > 1e-8) {
      stop(sprintf("%s scaling check FAILED: max diff = %.2e", label, max_diff))
    }
    message(sprintf("%s scaling check PASSED (max diff = %.2e)", label, max_diff))
  }

  check("Euclidean pairwise", proximity_euclidean_pairwise(X), "euclidean")
  check("Manhattan pairwise", proximity_manhattan_pairwise(X), "manhattan")

  invisible(TRUE)
}
