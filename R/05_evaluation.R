# =============================================================================
# 05_evaluation.R
#
# PURPOSE: Measure how well estimated cluster assignments match the true
#          cluster labels used to generate the data.
#
# Two metrics:
#
#   ARI — Adjusted Rand Index  (PRIMARY metric)
#     Compares two partitions by counting pairs of observations that are
#     placed in the SAME cluster or DIFFERENT clusters in both partitions.
#     The "adjusted" version corrects for the agreement expected by chance.
#     Range: −1 to 1.  Interpretation:
#       1.0  = perfect recovery of true clusters
#       0.0  = no better than random assignment
#      <0.0  = worse than random (unusual in practice)
#
#   Mean Silhouette Width  (SECONDARY metric)
#     Measures how well each observation fits its assigned cluster compared
#     to the next-best cluster, using the same distance matrix used for
#     clustering.  Range: −1 to 1.  Higher = cleaner, more compact clusters.
#     Unlike ARI, silhouette does NOT require knowing the true labels — it
#     evaluates the geometric quality of the clustering solution itself.
# =============================================================================


#' Adjusted Rand Index between true and estimated cluster assignments
#'
#' @param true_labels       Integer vector of true cluster memberships
#' @param estimated_labels  Integer vector of estimated cluster memberships
#' @return Scalar ARI value
compute_ari <- function(true_labels, estimated_labels) {
  mclust::adjustedRandIndex(true_labels, estimated_labels)
}


#' Mean silhouette width of a clustering solution
#'
#' @param estimated_labels  Integer vector of estimated cluster memberships
#' @param dist_matrix       The "dist" object used to produce the clustering
#' @return Scalar mean silhouette width, or NA if only 1 cluster was found
compute_silhouette <- function(estimated_labels, dist_matrix) {

  # Silhouette is undefined when all observations end up in one cluster
  n_unique_clusters <- length(unique(estimated_labels))
  if (n_unique_clusters < 2) return(NA_real_)

  sil <- cluster::silhouette(estimated_labels, dist_matrix)
  mean(sil[, "sil_width"])
}
