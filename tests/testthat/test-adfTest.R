library(testthat)
library(lumen)

# -------------------------------------------------------------------------
# Basic functionality
# -------------------------------------------------------------------------
test_that("adfTest: returns htest", {
  res <- suppressWarnings(adfTest(AirPassengers, lags = 3, type = "trend"))
  expect_s3_class(res, "htest")
})

test_that("adfTest: result has statistic, p.value and critical.values", {
  res <- suppressWarnings(adfTest(AirPassengers, lags = 3, type = "trend"))
  expect_false(is.null(res$statistic))
  expect_false(is.null(res$p.value))
  expect_false(is.null(res$critical.values))
})

test_that("adfTest: critical.values has 3 columns (1pct/5pct/10pct)", {
  res <- adfTest(AirPassengers, lags = 1, type = "drift")
  expect_equal(colnames(res$critical.values), c("1pct", "5pct", "10pct"))
})

test_that("adfTest: statistic count and names per type", {
  x <- diff(AirPassengers)

  r_none  <- suppressWarnings(adfTest(x, lags = 1, type = "none"))
  r_drift <- suppressWarnings(adfTest(x, lags = 1, type = "drift"))
  r_trend <- suppressWarnings(adfTest(AirPassengers, lags = 1, type = "trend"))

  expect_named(r_none$statistic,  "tau1")
  expect_named(r_drift$statistic, c("tau2", "phi1"))
  expect_named(r_trend$statistic, c("tau3", "phi2", "phi3"))
  expect_equal(rownames(r_trend$critical.values), c("tau3", "phi2", "phi3"))
})

# -------------------------------------------------------------------------
# Known reference values (identical to urca::ur.df, verified numerically)
# -------------------------------------------------------------------------
test_that("adfTest: AirPassengers reference values (trend, lags = 3)", {
  res <- suppressWarnings(adfTest(AirPassengers, lags = 3, type = "trend"))

  expect_equal(unname(res$statistic["tau3"]), -6.935821, tolerance = 1e-5)
  expect_equal(unname(res$statistic["phi2"]), 16.404216, tolerance = 1e-5)
  expect_equal(unname(res$statistic["phi3"]), 24.060097, tolerance = 1e-5)
  expect_equal(unname(res$parameter["lags"]), 3L)
})

test_that("adfTest: random walk reference values (drift)", {
  set.seed(1)
  rw <- cumsum(rnorm(100))

  res <- adfTest(rw, type = "drift")

  expect_equal(unname(res$statistic["tau2"]), -1.456120, tolerance = 1e-5)
  expect_equal(res$p.value, 0.516246, tolerance = 1e-5)
})

test_that("adfTest: matches tseries::adf.test (trend)", {
  skip_if_not_installed("tseries")

  set.seed(42)
  x <- cumsum(rnorm(150))
  for (k in c(0, 1, 3)) {
    res <- suppressWarnings(adfTest(x, lags = k, type = "trend"))
    ref <- suppressWarnings(tseries::adf.test(x, k = k))
    expect_equal(unname(res$statistic["tau3"]), unname(ref$statistic),
                 tolerance = 1e-8)
    expect_equal(res$p.value, ref$p.value, tolerance = 1e-8)
  }
})

# -------------------------------------------------------------------------
# Critical values (Dickey & Fuller 1981, Tables IV-VI; Fuller 1976)
# -------------------------------------------------------------------------
test_that("adfTest: critical values at n = 250 match the published tables", {
  # length 251 => 250 differences, hitting the n = 250 table row exactly;
  # phi3 guards against urca's transcription error (6.49/5.47 instead of
  # the published 6.34/5.39, DF 1981, Table VI)
  set.seed(7)
  y <- cumsum(rnorm(251))

  cv <- suppressWarnings(adfTest(y, type = "trend"))$critical.values

  expect_equal(unname(cv["tau3", ]), c(-3.99, -3.43, -3.13))
  expect_equal(unname(cv["phi2", ]), c(6.22, 4.75, 4.07))
  expect_equal(unname(cv["phi3", ]), c(8.43, 6.34, 5.39))
})

# -------------------------------------------------------------------------
# p-value behaviour
# -------------------------------------------------------------------------
test_that("adfTest: p-value in [0, 1]", {
  set.seed(1)
  res <- suppressWarnings(adfTest(cumsum(rnorm(80)), type = "drift"))
  expect_gte(res$p.value, 0)
  expect_lte(res$p.value, 1)
})

test_that("adfTest: stationary series rejects, with boundary warning", {
  set.seed(1)
  wn <- rnorm(200)

  expect_warning(
    res <- adfTest(wn, lags = 1, type = "drift"),
    "smaller than reported"
  )
  expect_lt(unname(res$statistic["tau2"]), -3)
  expect_equal(res$p.value, 0.01)
})

test_that("adfTest: random walk is not rejected", {
  set.seed(1)
  rw <- cumsum(rnorm(100))

  res <- adfTest(rw, lags = 1, type = "drift")

  expect_gt(res$p.value, 0.10)
})

# -------------------------------------------------------------------------
# Lag selection
# -------------------------------------------------------------------------
test_that("adfTest: parameter reports the lags actually used", {
  set.seed(3)
  x <- as.vector(arima.sim(list(ar = c(0.4, 0.3)), 200))

  r_fix <- suppressWarnings(adfTest(x, lags = 6, type = "drift"))
  r_aic <- suppressWarnings(adfTest(x, lags = 6, selectlags = "aic",
                                    type = "drift"))

  expect_equal(unname(r_fix$parameter["lags"]), 6L)
  expect_lte(unname(r_aic$parameter["lags"]), 6L)
})

test_that("adfTest: selectlags is case-insensitive", {
  set.seed(3)
  x <- as.vector(arima.sim(list(ar = c(0.4, 0.3)), 200))

  r_lower <- suppressWarnings(adfTest(x, lags = 6, selectlags = "aic"))
  r_upper <- suppressWarnings(adfTest(x, lags = 6, selectlags = "AIC"))

  expect_equal(r_lower$statistic, r_upper$statistic)
  expect_equal(r_lower$parameter, r_upper$parameter)
})

# -------------------------------------------------------------------------
# Interface
# -------------------------------------------------------------------------
test_that("adfTest: data.name and alternative are set", {
  set.seed(1)
  x <- cumsum(rnorm(50))

  r_drift <- adfTest(x, type = "drift")
  r_trend <- adfTest(x, type = "trend")

  expect_equal(r_drift$data.name, "x")
  expect_equal(r_drift$alternative, "stationary")
  expect_equal(r_trend$alternative, "trend-stationary")
})

test_that("adfTest: print.htest works", {
  set.seed(1)
  res <- adfTest(cumsum(rnorm(50)), type = "drift")
  expect_output(print(res), "Augmented Dickey-Fuller")
})

# -------------------------------------------------------------------------
# Input validation
# -------------------------------------------------------------------------
test_that("adfTest: NA in series throws error", {
  expect_error(adfTest(c(1, 2, NA, 4, 5), lags = 1))
})

test_that("adfTest: lags < 0 throws error", {
  expect_error(adfTest(rnorm(50), lags = -1))
})

test_that("adfTest: multivariate input throws error", {
  expect_error(adfTest(matrix(rnorm(40), ncol = 2)))
})

test_that("adfTest: too short series throws error", {
  expect_error(adfTest(rnorm(5), lags = 4))
})
