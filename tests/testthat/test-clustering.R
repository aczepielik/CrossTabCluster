test_that("hierarchical clustering works", {
  israeli_rows <- hclust_table(israeli_survey)
  israeli_cols <- hclust_table(israeli_survey, 2)

  expect_lte(
    abs(
      sum(israeli_rows$height) - sum(israeli_cols$height)
    ),
    1e-10
  )

  expect_lte(
    abs(
      sum(israeli_rows$height) - chisq.test(israeli_survey, simulate.p.value = TRUE)$statistic
    ),
    1e-10
  )

  ksarakil_rows <- hclust_table(ksarakil)
  ksarakil_cols <- hclust_table(ksarakil, 2)

  expect_lte(
    abs(
      sum(ksarakil_rows$height) - sum(ksarakil_cols$height)
    ),
    1e-10
  )

  expect_lte(
    abs(
      sum(ksarakil_rows$height) - chisq.test(ksarakil, simulate.p.value = TRUE)$statistic
    ),
    1e-10
  )

  expect_error(hclust_table(ksarakil, "row"))
  expect_error(hclust_table(ksarakil[1,], 1))
})
