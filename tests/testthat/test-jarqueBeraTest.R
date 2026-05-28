
# =============================================================================
# Tests for jarqueBeraTest()
# =============================================================================

library(testthat)

# --- Helpers ------------------------------------------------------------------

# Exact normal sample (large n, low Monte Carlo variance)
set.seed(42)
x_norm   <- rnorm(200)
x_skewed <- c(rexp(200))          # right-skewed
x_const  <- rep(1, 10)
x_na     <- c(rnorm(20), NA)
x_small  <- c(1, 2)               # n < 3

# =============================================================================
test_that("returns htest object with expected components", {
  out <- jarqueBeraTest(x_norm)
  expect_s3_class(out, "htest")
  expect_named(out$statistic, "X-squared")
  expect_true(is.numeric(out$p.value))
  expect_true(is.character(out$method))
  expect_true(is.character(out$data.name))
})

# --- p-value bounds -----------------------------------------------------------

test_that("p-value is in [0, 1] for chisq method", {
  out <- jarqueBeraTest(x_norm, method = "chisq")
  expect_gte(out$p.value, 0)
  expect_lte(out$p.value, 1)
})

test_that("p-value is in (0, 1] for mc method", {
  out <- jarqueBeraTest(x_norm, method = "mc", R = 499)
  expect_gt(out$p.value, 0)   # finite-sample correction ensures p > 0
  expect_lte(out$p.value, 1)
})

# --- Sensitivity --------------------------------------------------------------

test_that("clearly non-normal data yields small p-value (chisq)", {
  out <- jarqueBeraTest(x_skewed, robust = FALSE, method = "chisq")
  expect_lt(out$p.value, 0.05)
})

test_that("normal data yields large p-value (chisq)", {
  out <- jarqueBeraTest(x_norm, robust = FALSE, method = "chisq")
  expect_gt(out$p.value, 0.05)
})

# --- robust vs. classical -----------------------------------------------------

test_that("robust and classical return different statistics", {
  s_rob <- jarqueBeraTest(x_skewed, robust = TRUE,  method = "chisq")$statistic
  s_cls <- jarqueBeraTest(x_skewed, robust = FALSE, method = "chisq")$statistic
  expect_false(isTRUE(all.equal(s_rob, s_cls)))
})

# --- parameter field ----------------------------------------------------------

test_that("chisq method reports df = 2", {
  out <- jarqueBeraTest(x_norm, method = "chisq")
  expect_equal(out$parameter, c(df = 2))
})

test_that("mc method reports R in parameter", {
  out <- jarqueBeraTest(x_norm, method = "mc", R = 199)
  expect_equal(out$parameter, c(R = 199L))
})

# --- method name --------------------------------------------------------------

test_that("method string contains 'Robust' when robust = TRUE", {
  out <- jarqueBeraTest(x_norm, robust = TRUE)
  expect_match(out$method, "Robust")
})

test_that("method string does not contain 'Robust' when robust = FALSE", {
  out <- jarqueBeraTest(x_norm, robust = FALSE)
  expect_false(grepl("Robust", out$method))
})

test_that("method string reflects p-value method", {
  expect_match(jarqueBeraTest(x_norm, method = "chisq")$method, "Chi-square")
  expect_match(jarqueBeraTest(x_norm, method = "mc",    R = 99)$method,  "Monte Carlo")
})

# --- NA handling --------------------------------------------------------------

test_that("na.rm = TRUE removes NAs and runs", {
  expect_no_error(jarqueBeraTest(x_na, na.rm = TRUE))
})

test_that("na.rm = FALSE stops on NA", {
  expect_error(jarqueBeraTest(x_na, na.rm = FALSE), "missing values")
})

# --- Input validation ---------------------------------------------------------

test_that("constant vector throws error", {
  expect_error(jarqueBeraTest(x_const), "identical")
})

test_that("n < 3 throws error", {
  expect_error(jarqueBeraTest(x_small), "sample size")
})

test_that("matrix input throws error", {
  expect_error(jarqueBeraTest(matrix(1:9, 3, 3)), "numeric vector")
})

test_that("invalid R throws error", {
  expect_error(jarqueBeraTest(x_norm, method = "mc", R = 0),  "positive integer")
  expect_error(jarqueBeraTest(x_norm, method = "mc", R = -1), "positive integer")
})

test_that("data.name matches the variable name", {
  out <- jarqueBeraTest(x_norm)
  expect_equal(out$data.name, "x_norm")
})


