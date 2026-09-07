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

test_that("binomCIn: discrete methods are rejected", {
  # the root search needs non-integer x = p * n, which these methods
  # cannot deliver a smooth interval for
  for (m in c("mid-p", "blaker", "witting", "likelihood"))
    expect_error(binomCIn(p = 0.5, width = 0.1, method = m), "integer counts")
})

test_that("binomCIn: closed form methods can all be inverted", {
  # pratt is undefined at the lower end of the default interval (it needs
  # n * (1 - p) > 1), the search has to move the bracket up on its own
  for (m in c("wald", "wald-cc", "jeffreys", "jeffreys-mod",
              "clopper-pearson", "agresti-coull", "pratt", "arcsine",
              "logit", "wilson", "wilson-cc", "wilson-mod")) {
    n <- binomCIn(p = 0.4, width = 0.2, method = m)
    expect_length(n, 1L)
    expect_gt(n, 0)
  }
})

test_that("binomCIn: arguments are validated", {
  expect_error(binomCIn(p = 0, width = 0.1), "probability")
  expect_error(binomCIn(p = c(0.2, 0.3), width = 0.1), "probability")
  expect_error(binomCIn(p = 0.5, width = 0), "in \\(0, 1\\)")
  expect_error(binomCIn(p = 0.5, width = 0.1, conf.level = NA_real_),
               "must not be NA")
  expect_error(binomCIn(p = 0.5, width = 0.1, method = "nonesuch"))
  expect_error(binomCIn(p = 0.5, width = 0.1, interval = c(0, 10)), "interval")
})
