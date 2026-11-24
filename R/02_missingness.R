# =============================================================================
# 02_missingness.R
#
# PURPOSE: Introduce missing values into a complete dataset under two
#          missingness mechanisms:
#
#   MCAR — Missing Completely At Random
#          Each cell in the data matrix is independently set to NA with the
#          same probability (miss_rate). Missingness has no relationship to
#          any observed or unobserved value.
#
#   MAR  — Missing At Random
#          Whether a value is missing CAN be explained by observed information.
#          Here, observations with a HIGH value on a fully observed auxiliary
#          MAR driver have an elevated probability of missingness on the
#          clustering variables. The auxiliary driver is generated from a noisy
#          cluster-informative latent dimension, so MAR missingness can be
#          systematically heavier in some clusters than others. The driver is
#          not included in the clustering variables or proximity calculations.
#
# NOTE ON RATES: For both MCAR and MAR, miss_rate targets the overall
# proportion of missing cells among the clustering variables.
# =============================================================================


#' Top-level wrapper: introduce missing values by a chosen mechanism
#'
#' @param data       Complete numeric matrix (no NAs)
#' @param mechanism  "MCAR", "MAR", or "none" (returns data unchanged)
#' @param miss_rate  Proportion of values to set missing (0–1)
#' @param mar_driver Fully observed auxiliary variable used only for MAR
#'                   missingness. Required when mechanism = "MAR".
#'
#' @return The matrix with some values replaced by NA
introduce_missingness <- function(data, mechanism, miss_rate, mar_driver = NULL) {

  # Validation
  if (!mechanism %in% c("MCAR", "MAR", "none")) {
    stop("mechanism must be 'MCAR', 'MAR', or 'none'")
  }
  if (miss_rate < 0 || miss_rate > 1) {
    stop("miss_rate must be between 0 and 1")
  }

  # No missingness requested (e.g., complete-data baseline condition)
  if (mechanism == "none" || miss_rate == 0) return(data)

  if (mechanism == "MCAR") return(introduce_mcar(data, miss_rate))
  if (mechanism == "MAR")  return(introduce_mar(data, miss_rate, mar_driver))
}


# -----------------------------------------------------------------------------
# MCAR helper
# -----------------------------------------------------------------------------
# Every cell in the n × p matrix is independently set to NA with probability
# miss_rate. The missing pattern is unrelated to any data values.
introduce_mcar <- function(data, miss_rate) {
  data_out <- data

  for (j in seq_len(ncol(data_out))) {
    miss_mask <- runif(nrow(data_out)) < miss_rate
    data_out[miss_mask, j] <- NA
  }

  return(data_out)
}


# -----------------------------------------------------------------------------
# MAR helper
# -----------------------------------------------------------------------------
# A cell goes missing with probability that depends on the ROW's value of a
# fully observed auxiliary MAR driver. The driver is correlated with cluster
# membership but excluded from the analysis variables:
#
#   If driver > median(driver)  →  prob = miss_rate × 1.5  (higher risk)
#   If driver ≤ median(driver)  →  prob = miss_rate × 0.5  (lower risk)
#
# Because ~50% of rows are above the median, the AVERAGE probability across rows
# is approximately miss_rate. This is MAR because missingness depends on observed
# auxiliary information, not on the value that is made missing; because the
# driver is cluster-correlated, the missingness pattern can differ by true
# cluster.
introduce_mar <- function(data, miss_rate, mar_driver) {
  n <- nrow(data)
  data_out <- data

  if (is.null(mar_driver)) {
    stop("MAR missingness requires a fully observed mar_driver")
  }
  if (length(mar_driver) != n) {
    stop("mar_driver must have one value per row of data")
  }
  if (anyNA(mar_driver)) {
    stop("mar_driver must be fully observed")
  }

  # Per-row missingness probabilities for all clustering variables.
  high_driver <- mar_driver > median(mar_driver)
  rate_high <- min(miss_rate * 1.5, 0.95)  # cap at 95% to avoid all-missing rows
  rate_low  <- miss_rate * 0.5

  miss_probs <- ifelse(high_driver, rate_high, rate_low)  # length-n vector

  # Apply independently to every clustering variable.
  for (j in seq_len(ncol(data_out))) {
    miss_mask <- runif(n) < miss_probs
    data_out[miss_mask, j] <- NA
  }

  return(data_out)
}
