#' Calculate between profile psuedo distances
#'
#' @param profiles matrix with row profiles
#' @param D diagonal matrix with column weights
#' @param masses numeric vector with row masses
#'
#' @return square matrix with pseudo distances between profiles

calculate_psdist <- function(profiles, D, masses){
  k <- dim(profiles)[1]
  sp <- profiles %*% D %*% t(profiles) #scalar product
  distances <- matrix(rep(diag(sp), k), nrow = k, byrow = TRUE) +
    matrix(rep(diag(sp), k), nrow = k, byrow = FALSE) - 2*sp

  multipliers <- matrix(nrow = k, ncol = k)

  for (i in 1:k) {
    for (j in i:k) {
      multipliers[i,j] <- (masses[i] * masses[j])/(masses[i] + masses[j])
    }
  }

  multipliers[lower.tri(multipliers)] <- multipliers[upper.tri(multipliers)]
  distances*multipliers
}

#' Hierarchical clustering of contingency table
#'
#' Agglomerative clustring of rows or columns of contingency table in chi-square metric
#'
#' @param table contingency table to cluser. Must be coercible to matrix (e.g. class table)
#' @param dim chooses between clustering rows or columns. Character with dimension name (can be abbreviated)
#' or integer (1 for rows, 2 for columns)
#'
#' @return Objcect of class hclust. See ... for details
#'
#' @example
#' row_clust <- ct_hclust(israeli_worries)
#' plot(row_clust)
#'
#' @export

ct_hclust <- function(table, dim = "rows"){
  table <- as.matrix(table)
  N <- sum(table)
  freq_mat <- table/N
  row_masses <- rowSums(freq_mat)
  col_masses <- rowSums(t(freq_mat))

  if (grepl("r*", dim) | dim == 1) {
    profiles <- freq_mat/row_masses
    rownames(profiles) <- -seq_along(row_masses)
    D <- diag(1/col_masses)
    k <- length(row_masses)
    masses <- row_masses
  } else if (grepl("c*", dim) | dim == 2) {
    profiles <- t(freq_mat)/col_masses
    rownames(profiles) <- -seq_along(col_masses)
    D <- diag(1/row_masses)
    k <- length(col_masses)
    masses <- col_masses
  } else {
    stop("Neither rows nor colums selected to cluster \n Specify dim argument (see documentation)")
  }

  pseudo_dist <- calculate_psdist(profiles, D, masses)

  # Initialize clustering log------------
  merge <- matrix(nrow = k - 1, ncol = 2)
  height <- numeric(k - 1)

  # Actual clustering--------------------
  for (i in seq_along(height)) {
    labs <- rownames(profiles) #to be used in merging log

    to_cluster <- which(pseudo_dist == min(pseudo_dist[pseudo_dist > 0]),
                        arr.ind = TRUE) #choose pair with minimal pseudo distance
    to_cluster <- as.integer(to_cluster[1, ]) #select one of two symmetric results
    cluster_labs <- as.integer(labs[to_cluster])

    x <- to_cluster[1]
    y <- to_cluster[2]

    #Update log
    merge[i, ] <- cluster_labs
    height[i] <- pseudo_dist[x,y]

    #Construct new cluster profile
    cluster_mass <- masses[x] + masses[y]
    cluster_centroid <- (profiles[x, ]*masses[x] + profiles[y, ]*masses[y])/cluster_mass

    #Update profiles
    profiles <- profiles[-to_cluster, ]

    old_masses <- masses
    masses <- masses[-to_cluster]

    profiles <- rbind(profiles, cluster_centroid)

    rownames(profiles)[k - i] <- i
    #above doesn't work atthe last step where we have vector instead of matrix
    if (nrow(profiles) == 2) {
      rownames(profiles)[1] <- labs[-to_cluster]
    }

    masses <- c(masses, cluster_mass)

    #Update distance matrix
    if (i < k) {
      to_update <- setdiff(1:nrow(pseudo_dist), to_cluster)
      dist_to_cluster <- sapply(to_update,
                                function(ind)
                                  (old_masses[ind] + old_masses[x])/(old_masses[ind] + cluster_mass)*
                                  pseudo_dist[ind, x] +
                                  (old_masses[ind] + old_masses[y])/(old_masses[ind] + cluster_mass)*
                                  pseudo_dist[ind, y] -
                                  old_masses[ind]/(old_masses[ind] + cluster_mass)*
                                  pseudo_dist[x,y])
      pseudo_dist <- pseudo_dist[-to_cluster, -to_cluster]
      pseudo_dist <- cbind(rbind(pseudo_dist, dist_to_cluster), c(dist_to_cluster, 0))
      colnames(pseudo_dist) <- rownames(pseudo_dist) <- rownames(profiles)
    }
    #pseudo_dist <- calculate_psdist(profiles, D, masses)
  }

  #Building dendrogram-------------------
  indicators <- rev(as.integer(t(merge)))
  j <- 1
  while (j < length(indicators)) {
    if (indicators[j] < 0) {
      j <- j + 1
    } else{
      indicators <- append(indicators, merge[indicators[j], ], after = j)
      j <- j + 1
    }
  }
  order <- -unique(indicators[indicators < 0])

  result <- list(merge = merge, height = N*height, order = order, pseudo_dists = pseudo_dist)
  class(result) <- "hclust"

  result
}
