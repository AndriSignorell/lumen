library(testthat)
library(lumen)

# The p-value is interpolated from a table of critical values; outside that
# table kpssTest() warns that the true p-value is smaller/greater than the
# one reported (as tseries::kpss.test() does). Most tests here check other
# parts of the result, so they call it through this wrapper; the warning
# itself is tested explicitly at the end.
kpss <- function(...) suppressWarnings(kpssTest(...))

test_that("kpssTest: returns htest", {
  res <- kpss(AirPassengers, type = "mu")
  expect_s3_class(res, "htest")
})

test_that("kpssTest: statistic named KPSS is positive", {
  res <- kpss(AirPassengers, type = "mu")
  expect_named(res$statistic, "KPSS")
  expect_gt(unname(res$statistic), 0)
})

test_that("kpssTest: parameter named lags is non-negative integer", {
  res <- kpss(AirPassengers, type = "mu", lags = "short")
  expect_named(res$parameter, "lags")
  expect_true(res$parameter >= 0L)
  expect_true(is.integer(res$parameter))
})

test_that("kpssTest: critical.values matrix has correct shape", {
  res <- kpss(AirPassengers, type = "mu")
  cv  <- res$critical.values
  expect_true(is.matrix(cv))
  expect_equal(ncol(cv), 4L)
  expect_equal(colnames(cv), c("10pct", "5pct", "2.5pct", "1pct"))
})

test_that("kpssTest: type='mu' and type='tau' give different statistics", {
  res_mu  <- kpss(AirPassengers, type = "mu")
  res_tau <- kpss(AirPassengers, type = "tau")
  expect_false(isTRUE(all.equal(
    unname(res_mu$statistic),
    unname(res_tau$statistic)
  )))
})

test_that("kpssTest: lags='long' gives more lags than lags='short'", {
  res_s <- kpss(AirPassengers, lags = "short")
  res_l <- kpss(AirPassengers, lags = "long")
  expect_gte(unname(res_l$parameter), unname(res_s$parameter))
})

test_that("kpssTest: lags='nil' gives 0 lags", {
  res <- kpss(AirPassengers, lags = "nil")
  expect_equal(unname(res$parameter), 0L)
})

test_that("kpssTest: useLag overrides lags", {
  res <- kpss(AirPassengers, useLag = 3)
  expect_equal(unname(res$parameter), 3L)
})

test_that("kpssTest: stationary series gives small KPSS statistic", {
  set.seed(1)
  y   <- rnorm(200)           # iid ~ stationary
  res <- kpss(y, type = "mu", lags = "short")
  cv  <- res$critical.values[1, "10pct"]
  expect_lt(unname(res$statistic), cv)
})

test_that("kpssTest: random walk gives large KPSS statistic", {
  set.seed(2)
  y   <- cumsum(rnorm(2000))  # unit root - large n for reliable rejection
  res <- kpss(y, type = "mu", lags = "short")
  cv  <- res$critical.values[1, "1pct"]
  expect_gt(unname(res$statistic), cv)
})

test_that("kpssTest: NA values in input are silently removed", {
  y <- c(NA, as.numeric(AirPassengers), NA)
  expect_s3_class(kpss(y, type = "mu"), "htest")
})


test_that("kpssTest: warns when the p-value lies outside the table", {
  expect_warning(kpssTest(AirPassengers, type = "mu"),
                 "p-value smaller than reported")
  set.seed(1)
  expect_warning(kpssTest(rnorm(200), type = "mu", lags = "short"),
                 "p-value greater than reported")
})
