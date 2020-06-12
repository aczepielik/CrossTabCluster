sort_leaves <- function(merge, node) {
  if (merge[node, 1] < 0 & merge[node, 2] < 0) {
    c(merge[node, 1], merge[node, 2])
  } else if (merge[node, 1] < 0 & merge[node, 2] > 0) {
    c(sort_leaves(merge, merge[node, 2]), merge[node, 1])
  } else if (merge[node, 1] > 0 & merge[node, 2] < 0) {
    c(sort_leaves(merge, merge[node, 1]), merge[node, 2])
  } else if (merge[node, 1] > merge[node, 2]) {
    c(sort_leaves(merge, merge[node, 2]), sort_leaves(merge, merge[node, 1]))
  } else {
    c(sort_leaves(merge, merge[node, 1]), sort_leaves(merge, merge[node, 2]))
  }
}

lance_williams_update <- function(distances, masses, merge){
  inter_merge <- distances[merge[1], merge[2]]
  next_label <- max(0L, max(as.integer(rownames(distances)))) + 1L

  distances <- matrix(distances[, -merge], nrow = nrow(distances))

  old_distances <- matrix(distances[merge, ], nrow = 2)
  old_masses <- masses[merge]
  cluster_mass <- sum(old_masses)

  distances <- matrix(distances[-merge, ], ncol = ncol(distances))
  masses <- masses[-merge]

  new_distances <- sapply(seq_len(ncol(distances)), function(ind) {
    (masses[ind] + old_masses[1])/(masses[ind] + cluster_mass)*
      old_distances[1, ind] +
      (masses[ind] + old_masses[2])/(masses[ind] + cluster_mass)*
      old_distances[2, ind] -
      masses[ind]/(masses[ind] + cluster_mass)*
      inter_merge
  })

  new_distances <- cbind(rbind(distances, new_distances), c(new_distances, 0))
  colnames(new_distances)[ncol(new_distances)] <- next_label
  rownames(new_distances) <- colnames(new_distances)

  new_masses <- c(masses, cluster_mass)
  names(new_masses) <- colnames(new_distances)

  list(
    distances = new_distances,
    masses = new_masses
  )
}

#' Hierarchical clustering of contingency tables
#'
#' Applies Ward's agglomerative clustering to the rows or columns of the contingecy table.
#'
#' @inheritParams chisq_distance
#'
#' @return object of class hclust.
#' @export
#'
#' @examples
#' hclust_table(ksarakil)
#' hclust_table(israeli_survey, dimension = 2)
hclust_table <- function(table, dimension = 1) {
  if (dimension != 1 & dimension != 2) {
    stop("Wrong value for 'dimension' argument. Can be either 1 for rows or 2 for columns")
  }

  k <- dim(table)[dimension]

  if (k < 2) {
    stop("Must have two or more profiles to cluster")
  }

  #output prealocation
  merge <- matrix(nrow = k - 1, ncol = 2)
  height <- numeric(k - 1)

  #init
  distances <- cluster_distance(table, as.list(seq_len(k)), dimension)
  masses <- apply(table, dimension, sum)

  colnames(distances) <- rownames(distances) <- names(masses) <- -seq_len(k)

  for (i in seq_len(k-1)) {
    to_cluster <- as.integer(
      which(distances == min(distances[upper.tri(distances)]),
            arr.ind = TRUE)[1, ]
    )

    merge[i, ] <- as.integer(rownames(distances)[to_cluster])
    height[i] <- distances[to_cluster[1], to_cluster[2]]

    updates <- lance_williams_update(distances, masses, to_cluster)

    distances <- updates$distances

    masses <- updates$masses
  }

  structure(
    list(
      merge = merge,
      height = height,
      order = -sort_leaves(merge, nrow(merge)),
      labels = attr(table, "dimnames")[[dimension]],
      method = "ward.D2",
      call = match.call(),
      dist.method = "chisq"),
    class = "hclust")
}

simplified_ca <- function(table, dimension = 1) {
  if (dimension != 1 & dimension != 2) {
    stop("Wrong dimension specification. Can be either 1 for rows or 2 for columns")
  }

  if (dimension == 2) {
    table <- t(table)
  }

  N <- sum(table)

  row_masses <- rowSums(table)

  cc <- colSums(table)/N
  rr <- row_masses/N

  profiles <- table/row_masses
  decomposition <- svd(diag(sqrt(rr)) %*% sweep(profiles, 2, cc) %*% diag(sqrt(1/cc)))

  d <- decomposition$d
  n <- diag(1/sqrt(rr)) %*% decomposition$u
  m <- diag(sqrt(cc)) %*% decomposition$v
  f <- n %*% diag(d)
  g <- diag(d) %*% m

  list(f = f, n = n, g = g, n = n, m = m, d = d)
}

kmeans_table <- function(table, dimension = 1, k, ...) {
  args <- list(...)

  if ("x" %in% names(args)) {
    stop("Argument 'x' cannot be used as data is passed by kmeans_table itself")
  }

  if ("centers" %in% names(args)) {
    stop("Argument 'centers'  cannot be used as it is passed by kmeans_table itself")
  }

  if (dimension == 2) {
    table <- t(table)
  }

  row_masses <- apply(table, dimension, sum)
  indices <- sort(rep(seq_len(nrow(table)), row_masses))
  coordinates <- simplified_ca(table, dimension)$f
  coordinates <- coordinates[indices, 1:min(c(ncol(coordinates) - 1, k-1))]

  res <- do.call(kmeans, c(list(x = coordinates, centers = k), args))

  res$cluster.item <- res$cluster
  res$size.item <- res$size

  res$cluster <- as.integer(apply(table(indices, res$cluster), 1, which.max))
  res$size <- as.integer(table(res$cluster))

  res$totss.reduced <- res$totss
  res$withinss.reduced <- res$withinss
  res$betweenss.reduced <- res$betweenss

  res$withinss <- do.call(c, row_within_ss(table, clusters2list(res$cluster)))
  res$totss <- chisq.test(table, simulate.p.value = TRUE, B = 1)$statistic
  res$betweenss <- res$totss - sum(res$withinss)

  res
}
