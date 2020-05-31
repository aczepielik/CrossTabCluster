#' Upper bound for distance between clusters
#'
#'Under the null hypothesis of independence or homogeneity,
#'the maximal squared contrast in rows or columns of contingency table,
#'which can be interpreted as distance between rows is asympotically distributed as
#'the largest eigenvalue of the Wishart matrix W(I_min(I,J), max(I,J)), where I and J
#'are dimensions of the table.
#'
#' @param table object of class table. The contingency table to calculate the upper bound for.
#' @param alpha double. The returned value is 1-alpha quantile of the underlying distribution. Default is 0.05.
#' @param options list. Two elements taken into accout is "samples", which says how many repetitions to use
#' in simulating distribution (default: 1000) and "seed" used in sampling if not null (default: NULL).
#'
#' @return The 1-alpha quantile of the distribution of maximal estimated squared contrast
#' in tables of the size like the given table under the null hypothesis of independence or homogeneity.
#' @export
#'
#' @examples
#' distance_bound(israeli_survey, options = list(samples = 1000, seed = 42))
distance_bound <- function(table, alpha = 0.05, options = list(samples = 1000, seed = NULL)){
  a <- min(dim(table) - 1)
  b <- max(dim(table) - 1)

  if (!is.null(options[["seed"]])) {
    set.seed(options$seed)
  }

  stats::quantile(
    apply(stats::rWishart(options$samples, b, diag(a)), 3,
           function(mat) eigen(mat, only.values = TRUE)$values[1]),
    probs = 1 - alpha,
    names = FALSE)
}
