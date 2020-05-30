wishart_eigenvalue <- function(table, alpha, options = list(samples = 1000, seed = 42)){
  a <- min(dim(table) - 1)
  b <- max(dim(table) - 1)

  set.seed(options$seed)
  quantile(
    apply(rWishart(options$samples, b, diag(a)), 3,
           function(mat) eigen(mat, only.values = TRUE)$values[1]),
    probs = 1 - alpha,
    names = FALSE)
}
