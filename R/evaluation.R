sample_ind_mult_table <- function(n, row_probs, col_probs) {
  row_indices <- sample(seq_along(row_probs), size = n, replace = TRUE, prob = row_probs)
  col_indices <- sample(seq_along(col_probs), size = n, replace = TRUE, prob = col_probs)

  table(row_indices, col_indices)
}

sample_hom_prod_mult_table <- function(rowssums, col_probs) {
  M <- do.call(
    rbind,
    lapply(rowssums, function(n) {
      table(sample(seq_along(col_probs), size = n, replace = TRUE, prob = col_probs))
    })
  )

  rownames(M) <- seq_along(rowssums)
  as.table(M)
}

#' The gap statistic Work in progress
#'
#' Quantifies clustering quailty comparing it to siumulated clusters.
#' See Section 4.1 of the thesis for further details.
#' @inheritParams cluster_distance
#' @param func character. Name of clustering function. Currently "kmeans" and "hclust" are supported.
#' @param scheme character. Name of the null sampling scheme tables are genereated from. Can be either..
#' @param B integer. Number of sampled tables.
#'
#' @return list with elements: gap - value of the statistic and s -sclaled standard deviation
#'
#' @examples
#' clusters <- kmeans_table(israeli_survey, dimesnion = 1, k = 4)
#' gap_statistic(israeli_survey, dimension = 1, kmeans2list(clusters), "kmeans", "multinomial", 100)
gap_statistic <- function(table, dimension, clusters, func, scheme, B = 500) {
  if (dimension != 1 & dimension != 2) {
    stop("Wrong value for 'dimension' argument. Can be either 1 for rows or 2 for columns")
  }

  if (dimension == 2) {
    table <- t(table)
  }

  if(scheme == "multinomial") {
    samples <- lapply(seq_len(B),
                      function(i) {
                        sample_ind_mult_table(
                          sum(table), rowSums(table)/sum(table), colSums(table)/sum(table)
                        )})
  } else if (scheme == "product_multinomial" | scheme == "independent multinomial") {
    samples <- lapply(seq_len(B),
                      function(i) {
                        sample_hom_prod_mult_table(
                          rowSums(table), colSums(table)/sum(table)
                        )})
  } else {
    stop("Wrong scheme specification")
  }

  if (func == "kmeans" | func == "kmeans_table") {
    Wmc <- sapply(samples, function(tab) {
      sum(
        do.call(c,
                row_within_ss(
                  tab,
                  kmeans2list(kmeans_table(tab, 1, length(clusters)))
                )
        )
      )}
    )
  } else if (func == "hclust" | func == "hclust_table") {
    Wmc <- sapply(samples, function(tab) {
      sum(
        do.call(c,
        row_within_ss(
          tab,
          clusters2list(
            cutree(hclust_table(tab, 1), k = length(clusters))
          )
        )
      ))}
    )
  }

  W <- sum(do.call(c,row_within_ss(table, clusters)))
  list(gap = mean(log(Wmc)) - log(W), s = sqrt(1 + 1/B)*sd(log(Wmc)))
}

profile_distance_to_cluster <- function(distances, weights, profile_ind, cluster) {
  if (identical(profile_ind, cluster)) {
    0
  } else {
    dist_sum <- sum(distances[profile_ind, cluster] * weights[cluster])
    norm <-  ifelse(profile_ind %in% cluster,
                    sum(weights[cluster]) - weights[profile_ind],
                    sum(weights[cluster]))

    dist_sum/norm
  }
}

s <- function(a, b) {
  ifelse(a == 0 | a == b,
         0,
         ifelse(a < b, 1 - a/b, 1 - b/a))
}

#' Silhouette score
#'
#' Calculates silhouette scores for clustered profiles and mean score for the whole clustering result.
#'
#' @inheritParams chisq_distance
#' @return list with two elements. First element is a list with silhouette scores for each profile didvided
#' into clutsers, second one is a weighted mean of those scores.
#' @export
#'
#' @examples
#' clusters <- clusters2list(cutree(hclust_table(ksarakil, 1)))
#' silhouette(ksarakil, 1, clusters)
silhouette <- function(table, dimension, clusters) {
  if (dimension != 1 & dimension != 2) {
    stop("Wrong value for 'dimension' argument. Can be either 1 for rows or 2 for columns")
  }

  if (dimension == 2) {
    table <- t(table)
  }

  distances <- chisq_distance(table)
  rrn <- rowSums(table)

  a <- lapply(clusters, function(cl) {
    sapply(cl, function(profile_ind) {
      profile_distance_to_cluster(distances, rrn, profile_ind, cl)
    })
  })

  b <- lapply(clusters, function(cl) {
    sapply(cl, function(profile_ind) {
      min(
        sapply(clusters, function(cl_nest) {
          if (profile_ind %in% cl_nest) {
            NA
          } else {
            profile_distance_to_cluster(distances, rrn, profile_ind, cl_nest)
          }
        }),
        na.rm = TRUE)
    })
  })

  silhouettes <-  lapply(seq_along(clusters), function(i) s(a[[i]], b[[i]]))

  list(
    silhouettes = silhouettes,
    score = sum(do.call(c, silhouettes) * rrn[do.call(c, clusters)])/sum(rrn)
  )
}

aic_multinomial <- function(table, dimension, clusters) {
  if (dimension == 2) {
    table <- t(table)
  }

  const <- lgamma(sum(table) + 1) -
    sum(lgamma(table + 1)) +
    sum(rowSums(table) * log(rowSums(table)/sum(table)))

  cls <- lapply(clusters, function(cl) matrix(slice_table(table, cl, 1), ncol = ncol(table)))

  -2* (const + sum(sapply(cls, function(cl) sum(colSums(cl) * log(colSums(cl)/sum(cl)), na.rm = TRUE) )) ) +
    2*((ncol(table) - 1)*length(clusters) + nrow(clusters) - 1)
}

aic_prod_mult <- function(table, dimension, clusters) {
  if (dimension == 2) {
    table <- t(table)
  }

  const <- sum(lgamma(rowSums(table) + 1)) -
    sum(lgamma(table + 1))

  cls <- lapply(clusters, function(cl) matrix(slice_table(table, cl, 1), ncol = ncol(table)))

  -2* (const + sum(sapply(cls, function(cl) sum(colSums(cl) * log(colSums(cl)/sum(cl)), na.rm=TRUE) ) ) ) +
    2*(ncol(table) - 1)*length(clusters)
}

#' Akaike Information Criterion
#'
#' See section 4.1 in thesis, for details
#' @inheritParams cluster_distance
#' @param scheme character. Name of the assumed table's generation scheme.
#' Can be either "multinomial" or "independent multinomial"
#'
#' @return double. Value of AIC statisctic.
#' @export
#'
#' @examples
#' clusters <- clusters2list(cutree(hclust_table(ksarakil, 1)))
#' aic(ksarakil, 1, clusters, "independent multinomial")
aic <- function(table, dimension, clusters, scheme) {
  if (dimension != 1 & dimension != 2) {
    stop("Wrong value for 'dimension' argument. Can be either 1 for rows or 2 for columns")
  }

  if(scheme == "multinomial") {
    aic_multinomial(table, dimension, clusters)
  } else if (scheme == "product_multinomial" | scheme == "independent multinomial") {
    aic_prod_mult(table, dimension, clusters)
  } else {
    stop("Wrong scheme specification")
  }
}
