test_that("killconf works properly and show error if needed", {
  expect_error({
    data(lalonde, package = "Matching")

    data <- dplyr::mutate(lalonde,
      black = factor(black),
      hisp = factor(hisp),
      married = factor(married),
      treat = as.character(treat)
    )

    killconf(data = data, treatment = "treat", indvar = c("age", "educ", "black", "hisp", "married", "re75", "re74"), outvar = "re78", categorical = FALSE, R = 5)
  })
  expect_error({
    data(lalonde, package = "Matching")

    data <- dplyr::mutate(lalonde,
      black = factor(black),
      hisp = factor(hisp),
      married = factor(married)
    )

    killconf(data = data, treatment = "treat", indvar = c("age", "educ", "black", "hisp", "married", "re75", "re74"), outvar = "re78", categorical = TRUE, R = 5)
  })
})
