test_that("cluster correcteness works", {
  expect_warning(
    test1 <- check_cluster_correctness(ksarakil, list(1, c(1, 2), c(3, 4, 5)))
  )

  expect_false(test1)

  expect_false(
    check_cluster_correctness(ksarakil, list(3, c(4, 5), c(6, 7, 8)), dimension = 2)
    )

  expect_true(
    check_cluster_correctness(ksarakil, list(3, c(4, 5), c(6, 7, 8)), dimension = 1)
  )
})

test_that("slicing", {
  expect_equal(slice_table(ksarakil, c(1:3), 1), ksarakil[1:3, ])
  expect_equal(slice_table(ksarakil, c(1:3), 2), t(ksarakil[, 1:3]))

  expect_error(slice_table(ksarakil, c(1:3), 'row'))
})

test_that("Gramm Schmidt", {
  set.seed(42)

  M <- matrix(runif(20), nrow = 5, ncol = 4)
  gsm <- gramm_schmidt(M)
  expect_true(all(t(gsm)%*%gsm - diag(4) < 1e-10))

  gsm2 <- gramm_schmidt(M, normalize = FALSE)
  expect_true(all(t(gsm)%*%gsm - diag(t(gsm)%*%gsm) < 1e-10))
})

test_that("unit norm", {
  set.seed(42)
  x <- runif(5)
  expect_false(check_unit_norm(x))
  expect_true(check_unit_norm(x/sqrt(sum(x^2))))
})

test_that("orthogonal complement", {
  set.seed(42)
  x <- runif(5)
  expect_error(orthonormal_complement(x))
  expect_equal(dim(orthonormal_complement(x/sqrt(sum(x^2)))), c(5, 4))
  expect_true(all(apply(orthonormal_complement(x/sqrt(sum(x^2))), 2, check_unit_norm)))
  expect_true(all(apply(orthonormal_complement(x/sqrt(sum(x^2))), 2,
                        function(vec) t(x/sqrt(sum(x^2))) %*% vec <= 1e-10)))
})

test_that("l2 marginals work", {
  expect_equal(row_sums_l2(ksarakil), col_sums_l2(t(ksarakil)))
})

test_that("chi square distances sum up to chi-square statistic", {
  colSums(israeli_survey) -> centroid
  distances <- chisq_distance(rbind(israeli_survey, centroid))

  expect_lte(
    abs(
      sum(rowSums(israeli_survey) * distances[9, 1:8]) -
        chisq.test(israeli_survey, simulate.p.value = TRUE)$statistic
    ),
    1e-10
  )

  expect_error(chisq_distance(israeli_survey, as.list(1:8), 2))
})

test_that("cluster distances are like in hirotsu 2009", {
  dists <- cluster_distance(israeli_survey, as.list(seq_len(nrow(israeli_survey))), 1)
  expect_lte(abs(dists[1, 2] - 13.79), 0.01)
  expect_lte(abs(dists[7, 5] - 2.06), 0.01)
  expect_lte(abs(dists[5, 6] - 0.13), 0.01)

  expect_error(cluster_distance(israeli_survey, as.list(1:8), 2))
})

test_that("generalized cluster distances are like in hirotsu 2009", {
  expect_lte(
    abs(
      generalized_distance(israeli_survey, list(c(1, 2, 3, 5, 4, 6, 7), 8), 1) -
        77.9),
    0.01)

  expect_lte(
    abs(
      generalized_distance(israeli_survey, list(c(2, 3, 5, 4, 6, 7), 1), 1) -
        20.77),
    0.01)

  expect_lte(
    abs(
      generalized_distance(israeli_survey, list(c(2, 4), c(5, 3, 6, 7), 1), 1) -
        21.76),
    0.01)

  expect_lte(
    abs(
      generalized_distance(israeli_survey, list(c(7, 4, 3), c(5, 6), 1, 2), 1) -
        25.36),
    0.01)

  expect_equal(
    generalized_distance(israeli_survey, as.list(1:8), 1),
    generalized_distance(israeli_survey, as.list(1:5), 2)
    )

  expect_error(generalized_distance(israeli_survey, as.list(1:8), 2))
})
