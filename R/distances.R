check_cluster_correctness <- function(table, clusters, dimension = 1) {
  size <- dim(table)[dimension]
  flattened <- do.call(c, clusters)

  all(unique(flattened) == flattened) & max(flattened) <= size
}

slice_table <- function(table, indices, dimension) {
  if (dimension == 1) {
    profiles <- table[indices, ]
  } else if (dimension == 2) {
    profiles <- t(table[, indices])
  } else {
    stop("Wrong value for 'dimension' argument. Can be either 1 for rows or 2 for columns")
  }
  profiles
}

ward_coeff_xy <- function(x, y) {
  x*y/(x+y)
}

ward_coeffs <- function(x) {
  tmp <- expand.grid(x, x)
  matrix(ward_coeff_xy(tmp$Var1, tmp$Var2), nrow = length(x))
}

proj <- function(u, v) {
  c((t(v) %*% u)/(t(u) %*% u))* u
}

gramm_schmidt <- function(M, normalize = TRUE) {
  U <- M

  for (col_ind in 2L:ncol(U)) {
      U[, col_ind] <- U[, col_ind] -
        rowSums(sapply(1L:(col_ind -1), function(i) proj(U[, i], U[, col_ind])))
  }

  if (normalize) {
    apply(U, 2, function(col) col/sqrt(sum(col^2)))
  } else {
    U
  }
}

check_unit_norm <- function(vector, tol = 1e-05){
  abs(sum(vector^2) - 1) < tol
}

orthonormal_complement <- function(vector) {
  stopifnot(check_unit_norm(vector))

  base <- diag(length(vector))
  base[,1] <- vector

  gramm_schmidt(base)[, -1]
}

row_sums_l2 <- function(table) {
  sqrt(rowSums(table))/sqrt(sum(table))
}

col_sums_l2 <- function(table) {
  sqrt(colSums(table))/sqrt(sum(table))
}

row_distance_W <- function(table, clusters) {
  sapply(clusters, function(clt) {
    1/sqrt((sum(table[clt, ])/sum(table))) *
      t(orthonormal_complement(col_sums_l2(table))) %*%
      (
        colSums(matrix(table[clt, ], ncol = ncol(table))) *
          1/sqrt(colSums(table))
      )
  })
}

#' Matrix of chi-square distances
#'
#' Calculates chi-square distances between selected rows or columns of the contingency table.
#'
#' @param table object of class "table".
#' @param indices integer vector. Which rows or columns should the function calculate distances between.\
#' Default is all the row or all the columns depending on what is chosen in the argument "dimension".
#' @param dimension integer. Whether to calculate distances between rows (1) or columns (2). Default is 1.
#'
#' @return Matrix of size length(indices) x length(indices) containing chi-square distances between selected
#' rows or columns.
#'
#' @export
#'
#' @examples
#' data(ksarakil)
#' chisq_distance(ksarakil)
chisq_distance <- function(table, indices = seq_len(dim(table)[dimension]), dimension = 1) {
  profiles <- slice_table(table, indices, dimension)/apply(table, dimension, sum)

  coordinate_weights <- apply(table, ifelse(dimension == 1, 2, 1), sum)/sum(table)

  D <- diag(1/coordinate_weights)
  k <- nrow(profiles)

  sp <- profiles %*% D %*% t(profiles) #scalar product

  matrix(rep(diag(sp), k), nrow = k, byrow = TRUE) +
    matrix(rep(diag(sp), k), nrow = k, byrow = FALSE) - 2*sp
}

#' Matrix of Ward's distances between clusters
#'
#' Calculates chi-square distances between selected rows or columns of the contingency table.
#'
#' @inheritParams chisq_distance
#' @param clusters list of integer vectors. Each vector should define a cluster
#' by specifing row or column indices of its memebrs. Clusters must not overalap.
#'
#' @return Matrix of size length(clusters) x length(clusters) containing distances between selected clusters of
#' rows or columns.
#'
#' @export
#'
#' @examples
#' data(israeli_survey)
#' cluster_distance(israeli_survey, as.list(seq_len(nrow(israeli_survey))), 1)
#'
#' cluster_distance(israeli_survey, list(1, 2, c(3, 5), c(4, 6, 7), 8), 1)
cluster_distance <- function(table, clusters, dimension = 1) {
  if (!check_cluster_correctness(table, clusters, dimension)) {
    stop("Wrong cluster specification. Clusters overlap or indices or outside of table range")
  }

  rest_of_profiles <- setdiff(seq_len(dim(table)[dimension]), do.call(c, clusters))

  collapsed_table <- do.call(
    rbind,
    lapply(
      c(clusters, rest_of_profiles),
      function(cl) colSums(matrix(slice_table(table, cl, dimension), ncol = dim(table)[3 - dimension]))
    )
  )

  rest_cluster <- ifelse(length(rest_of_profiles) == 0, 0, 1)
  row_masses <- rowSums(collapsed_table)[seq_len(nrow(collapsed_table) - rest_cluster)]


  chisq_distance(collapsed_table)[seq_len(nrow(collapsed_table) - rest_cluster),
                                  seq_len(nrow(collapsed_table) - rest_cluster)] *
    ward_coeffs(row_masses)
}

#' Generalized distance between clusters
#'
#' @inheritParams cluster_distance
#'
#' @return double. Value of generalized distance between clusters
#' @export
#'
#' @examples
#' data(israeli_survey)
#' generalized_distance(israeli_survey, as.list(seq_len(nrow(israeli_survey))), 1)
#'
#' generalized_distance(israeli_survey, list(1, 2, c(3, 5), c(4, 6, 7)), 1)
generalized_distance <- function(table, clusters, dimension = 1) {
  if (!check_cluster_correctness(table, clusters, dimension)) {
    stop("Wrong cluster specification. Clusters overlap or indices or outside of table range")
  }

  if(dimension == 2) {
    table <- t(table)
  }

  W <- row_distance_W(table, clusters)
  Tk <- sapply(clusters, function(clt) sqrt(sum(table[clt, ])))
  Tk_normed <- Tk/sqrt(sum(Tk^2))
  Tk_comp <- orthonormal_complement(Tk_normed)

  eigen(t(Tk_comp) %*% t(W) %*% W %*% Tk_comp, only.values = TRUE)$values[1]
}
