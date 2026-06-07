library(testthat)
library(lumen)

test_that("meanCIn: returns positive numeric", {
  n <- meanCIn(ci = c(4, 6), sd = 2)
  expect_type(n, "double")
  expect_gt(n, 0)
})

test_that("meanCIn: achieved CI width matches target", {
  target_half <- 1   # ci = c(4,6) -> half-width = 1
  n <- ceiling(meanCIn(ci = c(4, 6), sd = 2))
  # t-based CI half-width: qt(0.975, n-1) * sd / sqrt(n)
  achieved <- qt(0.975, n - 1) * 2 / sqrt(n)
  expect_equal(achieved, target_half, tolerance = 0.05)
})

test_that("meanCIn: larger sd requires larger n", {
  n1 <- meanCIn(ci = c(4, 6), sd = 1)
  n2 <- meanCIn(ci = c(4, 6), sd = 2)
  expect_gt(n2, n1)
})

test_that("meanCIn: narrower CI requires larger n", {
  n_narrow <- meanCIn(ci = c(4.5, 5.5), sd = 2)
  n_wide   <- meanCIn(ci = c(4,   6  ), sd = 2)
  expect_gt(n_narrow, n_wide)
})

test_that("meanCIn: higher conf.level requires larger n", {
  n95 <- meanCIn(ci = c(4, 6), sd = 2, conf.level = 0.95)
  n99 <- meanCIn(ci = c(4, 6), sd = 2, conf.level = 0.99)
  expect_gt(n99, n95)
})

test_that("meanCIn: norm=TRUE uses z instead of t (gives smaller n)", {
  n_t <- meanCIn(ci = c(4, 6), sd = 2, norm = FALSE)
  n_z <- meanCIn(ci = c(4, 6), sd = 2, norm = TRUE)
  expect_gt(n_t, n_z)
})
