check_cluster_correctness <- function(table, clusters, dimension = 1) {
  size <- dim(table)[dimension]
  flattened <- do.call(c, clusters)

  unique(flattened) == flattened & max(flattened) <= size
}

chisq_distance <- function(table, indices, dimension = 1) {
  profiles <- table[dimension, indices]/apply(table, dimension, sum)
  weights <- apply(table, ifelse(dimension == 1, 2, 1), sum)
  masses <- apply(table, dimension, sum)

  profiles
}

cluster_distance <- function(table, clusters, dimension = 1) {
  NULL
}

generealized_distance <- function(table, clusters, dimension = 1) {
  NULL
}
