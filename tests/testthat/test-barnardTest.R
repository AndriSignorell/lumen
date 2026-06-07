library(testthat)
library(lumen)

# Mehta et al. (2003) example
tab_mehta <- matrix(c(7, 12, 8, 3), nrow = 2,
                    dimnames = list(treat = c("vaccine","placebo"),
                                    infection = c("yes","no")))

# Small balanced table
tab_small <- matrix(c(8, 14, 1, 3), nrow = 2)

test_that("barnardTest: returns htest", {
  res <- barnardTest(tab_small)
  expect_s3_class(res, "htest")
})

test_that("barnardTest: p.value in [0,1]", {
  res <- barnardTest(tab_small)
  expect_true(res$p.value >= 0 && res$p.value <= 1)
})

test_that("barnardTest: two-sided p >= one-sided p", {
  res_two  <- barnardTest(tab_mehta, alternative = "two.sided")
  res_less <- barnardTest(tab_mehta, alternative = "less")
  expect_gte(res_two$p.value, res_less$p.value - 1e-6)
})

test_that("barnardTest: non-2x2 matrix throws error", {
  tab3x2 <- matrix(1:6, nrow = 3)
  expect_error(barnardTest(tab3x2))
})

test_that("barnardTest: Mehta example one-sided p < 0.05", {
  res <- barnardTest(tab_mehta, alternative = "less")
  expect_lt(res$p.value, 0.05)
})

test_that("barnardTest: result has method field", {
  res <- barnardTest(tab_small)
  expect_false(is.null(res$method))
})

test_that("barnardTest: fixed=2 (column margins) works", {
  res <- barnardTest(tab_small, fixed = 2)
  expect_s3_class(res, "htest")
})
