library(testthat)
library(lumen)

test_that("adfTest: returns htest", {
  res <- adfTest(AirPassengers, lags = 3, type = "trend")
  expect_s3_class(res, "htest")
})

test_that("adfTest: result has statistic and critical.values", {
  res <- adfTest(AirPassengers, lags = 3, type = "trend")
  expect_false(is.null(res$statistic))
  expect_false(is.null(res$critical.values))
})

test_that("adfTest: critical.values has 3 columns (1pct/5pct/10pct)", {
  res <- adfTest(AirPassengers, lags = 1, type = "drift")
  expect_equal(colnames(res$critical.values), c("1pct", "5pct", "10pct"))
})

test_that("adfTest: type='none' gives 1 statistic value", {
  res <- adfTest(diff(AirPassengers), lags = 1, type = "none")
  expect_equal(length(res$statistic), 1L)
})

test_that("adfTest: type='drift' gives 2 statistic values", {
  res <- adfTest(diff(AirPassengers), lags = 1, type = "drift")
  expect_equal(length(res$statistic), 2L)
})

test_that("adfTest: type='trend' gives 3 statistic values", {
  res <- adfTest(AirPassengers, lags = 1, type = "trend")
  expect_equal(length(res$statistic), 3L)
})

test_that("adfTest: stationary series gives strongly negative tau", {
  # white noise is stationary; ADF tau should be very negative
  set.seed(1)
  wn <- rnorm(200)
  res <- adfTest(wn, lags = 1, type = "drift")
  tau <- res$statistic[1]
  expect_lt(tau, -3)
})

test_that("adfTest: random walk gives tau near critical value or above", {
  set.seed(1)
  rw <- cumsum(rnorm(100))
  res <- adfTest(rw, lags = 1, type = "drift")
  # tau for a unit root should typically be > -3.5 (not strongly reject)
  tau <- res$statistic[1]
  expect_gt(tau, -6)  # loose bound: just checks it's a plausible unit-root value
})

test_that("adfTest: NA in series throws error", {
  expect_error(adfTest(c(1, 2, NA, 4, 5), lags = 1))
})

test_that("adfTest: lags < 0 throws error", {
  expect_error(adfTest(rnorm(50), lags = -1))
})
