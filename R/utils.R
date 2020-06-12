#' Refactor cluster memebership
#'
#' Turn clsuter membershiip vector either from kmeans or cutree into list of vectors of cluster members
#'
#' @param vec vector of cluster membership indices.
#' @param kmeans_result object of class kmeans
#'
#' @return list of vectors of cluster members
#' @export
#'
#' @examples
#' clusters2list(c(1, 1, 2, 3, 1, 2, 3))
#'
#' clusters <- kmeans_table(israeli_survey, 1, k = 4)
#' kmeans2list(clusters)
clusters2list <- function(vec) {
  res <- as.list(unname(unclass(by(seq_along(vec), vec, c))))
  attr(res, "call") <- NULL
  res
}

#' @rdname clusters2list
#' @export
kmeans2list <- function(kmeans_result) {
  clusters2list(kmeans_result$cluster)
}

#' Within-cluster sum of squares
#'
#' Calclulates sum of squared distances between row profiles and their cluster centroids
#' @param table
#' @param clusters
#'
#' @return
#' @export
#'
#' @examples
row_within_ss <- function(table, clusters) {
  n <- sum(table)
  cc <- colSums(table)/n
  rrn <- rowSums(table)
  #rr <- rrn/n

  profiles <- sweep(table, 1, rrn, `/`)
  lapply(clusters, function(cl) {
    subprofiles <- matrix(profiles[cl, ], ncol = ncol(table))
    centroid <- colSums(
      matrix(sweep(subprofiles, 1, rrn[cl], `*`), ncol = ncol(table)))/sum(rrn[cl]
      )
    sum(rrn[cl] *
          rowSums(
            t(
              apply(sweep(subprofiles, 2, centroid), 1, function(profile) profile/sqrt(cc)))^2))
  })

}
