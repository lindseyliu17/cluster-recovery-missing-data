# =============================================================================
# 04_clustering.R
#
# PURPOSE: Apply clustering algorithms to a precomputed distance matrix
#          and return a vector of cluster assignments.
#
# Both methods work directly from a dissimilarity matrix, which is important
# because some proximity matrices (e.g., pairwise-complete) may not correspond
# to any Euclidean geometry — you cannot always recover coordinates from them.
#
# Methods implemented:
#
#   PAM (Partitioning Around Medoids)
#     Non-hierarchical. Selects k actual observations as "medoids" (cluster
#     centers) and assigns every other observation to its nearest medoid.
#     Unlike k-means, PAM works with arbitrary distances (not just Euclidean)
#     and is more robust to outliers.  It accepts precomputed distance matrices
#     directly, making it ideal for this simulation.
#
#   Average Linkage Hierarchical Clustering (UPGMA)
#     Hierarchical. At each step it merges the two clusters whose average
#     pairwise distance is smallest, weighting the update by the number of
#     observations in each cluster. This is the observation-weighted average
#     linkage version, not McQuitty/WPGMA's equal-cluster-weight update.
#     In R, UPGMA is implemented as hclust(method = "average").
# =============================================================================


#' Apply a clustering method to a precomputed distance matrix
#'
#' @param dist_matrix  A "dist" object (from stats::dist, cluster::daisy, etc.)
#' @param method       "pam" or "average"
#' @param k            Number of clusters to extract (default = 3)
#'
#' @return Integer vector of cluster assignments (length = n observations)
cluster_from_distance <- function(dist_matrix, method, k = 3) {

  if (!method %in% c("pam", "average")) {
    stop("method must be 'pam' or 'average'")
  }

  if (method == "pam")     return(cluster_pam(dist_matrix, k))
  if (method == "average") return(cluster_average_linkage(dist_matrix, k))
}


# -----------------------------------------------------------------------------
# PAM
# -----------------------------------------------------------------------------
#' Partitioning Around Medoids
#'
#' @param dist_matrix  "dist" object
#' @param k            Number of clusters
#' @return Integer vector of cluster assignments
cluster_pam <- function(dist_matrix, k) {
  fit <- cluster::pam(x = dist_matrix, k = k, diss = TRUE)
  as.integer(fit$clustering)
}


# -----------------------------------------------------------------------------
# Average Linkage Hierarchical Clustering
# -----------------------------------------------------------------------------
#' Average linkage (UPGMA) hierarchical clustering
#'
#' Builds a complete dendrogram with hclust(method = "average"), then cuts it
#' into k groups. This is UPGMA, which weights linkage updates by cluster size.
#'
#' @param dist_matrix  "dist" object
#' @param k            Number of clusters to cut the dendrogram into
#' @return Integer vector of cluster assignments
cluster_average_linkage <- function(dist_matrix, k) {
  hc <- hclust(dist_matrix, method = "average")
  as.integer(cutree(hc, k = k))
}
