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
      order = seq_len(k),
      labels = attr(table, "dimnames")[[dimension]],
      method = "ward.D2",
      call = match.call(),
      dist.method = "chisq"),
    class = "hclust")
}

kmeans_table <- function(table, dimension) {
  NULL
}
