test_that("upper bound works like in hirotsu", {
  expect_lte(
    abs(
      distance_bound(israeli_survey, options = list(samples = 1000, seed = 42)) - 23.55
    ),
    0.1
  )
})
