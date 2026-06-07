library(testthat)
library(lumen)

test_that("binomCIn: returns a single numeric", {
  n <- binomCIn(p = 0.5, width = 0.1)
  expect_type(n, "double")
  expect_length(n, 1L)
})

test_that("binomCIn: achieved CI width ~ target width", {
  target <- 0.1
  n <- ceiling(binomCIn(p = 0.5, width = target))
  actual_width <- unname(diff(binomCI(x = round(0.5 * n), n = n)[-1]))
  expect_equal(actual_width, target, tolerance = 0.01)
})

test_that("binomCIn: larger width requires smaller n", {
  n_narrow <- binomCIn(p = 0.5, width = 0.05)
  n_wide   <- binomCIn(p = 0.5, width = 0.10)
  expect_gt(n_narrow, n_wide)
})

test_that("binomCIn: p=0.5 is worst case (largest n)", {
  n_half    <- binomCIn(p = 0.5, width = 0.1)
  n_extreme <- binomCIn(p = 0.1, width = 0.1)
  expect_gt(n_half, n_extreme)
})

test_that("binomCIn: higher conf.level requires larger n", {
  n95 <- binomCIn(p = 0.5, width = 0.1, conf.level = 0.95)
  n99 <- binomCIn(p = 0.5, width = 0.1, conf.level = 0.99)
  expect_gt(n99, n95)
})

test_that("binomCIn: result is positive", {
  expect_gt(binomCIn(p = 0.3, width = 0.08), 0)
})
