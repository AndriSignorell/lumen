library(testthat)
library(lumen)

test_that("kpssTest: returns htest", {
  res <- kpssTest(AirPassengers, type = "mu")
  expect_s3_class(res, "htest")
})

test_that("kpssTest: statistic named KPSS is positive", {
  res <- kpssTest(AirPassengers, type = "mu")
  expect_named(res$statistic, "KPSS")
  expect_gt(unname(res$statistic), 0)
})

test_that("kpssTest: parameter named lags is non-negative integer", {
  res <- kpssTest(AirPassengers, type = "mu", lags = "short")
  expect_named(res$parameter, "lags")
  expect_true(res$parameter >= 0L)
  expect_true(is.integer(res$parameter))
})

test_that("kpssTest: critical.values matrix has correct shape", {
  res <- kpssTest(AirPassengers, type = "mu")
  cv  <- res$critical.values
  expect_true(is.matrix(cv))
  expect_equal(ncol(cv), 4L)
  expect_equal(colnames(cv), c("10pct", "5pct", "2.5pct", "1pct"))
})

test_that("kpssTest: type='mu' and type='tau' give different statistics", {
  res_mu  <- kpssTest(AirPassengers, type = "mu")
  res_tau <- kpssTest(AirPassengers, type = "tau")
  expect_false(isTRUE(all.equal(
    unname(res_mu$statistic),
    unname(res_tau$statistic)
  )))
})

test_that("kpssTest: lags='long' gives more lags than lags='short'", {
  res_s <- kpssTest(AirPassengers, lags = "short")
  res_l <- kpssTest(AirPassengers, lags = "long")
  expect_gte(unname(res_l$parameter), unname(res_s$parameter))
})

test_that("kpssTest: lags='nil' gives 0 lags", {
  res <- kpssTest(AirPassengers, lags = "nil")
  expect_equal(unname(res$parameter), 0L)
})

test_that("kpssTest: useLag overrides lags", {
  res <- kpssTest(AirPassengers, useLag = 3)
  expect_equal(unname(res$parameter), 3L)
})

test_that("kpssTest: stationary series gives small KPSS statistic", {
  set.seed(1)
  y   <- rnorm(200)           # iid ~ stationary
  res <- kpssTest(y, type = "mu", lags = "short")
  cv  <- res$critical.values[1, "10pct"]
  expect_lt(unname(res$statistic), cv)
})

test_that("kpssTest: random walk gives large KPSS statistic", {
  set.seed(2)
  y   <- cumsum(rnorm(2000))  # unit root - large n for reliable rejection
  res <- kpssTest(y, type = "mu", lags = "short")
  cv  <- res$critical.values[1, "1pct"]
  expect_gt(unname(res$statistic), cv)
})

test_that("kpssTest: NA values in input are silently removed", {
  y <- c(NA, as.numeric(AirPassengers), NA)
  expect_s3_class(kpssTest(y, type = "mu"), "htest")
})
