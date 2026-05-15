
# -------------------------------------------------------------------------
# Basic functionality
# -------------------------------------------------------------------------
test_that("bartelsRankTest returns an htest object", {
  
  set.seed(1)
  x <- rnorm(20)
  
  res <- bartelsRankTest(x)
  
  expect_s3_class(res, "htest")
  expect_named(res$statistic, "z")
  expect_named(res$parameter, "n")
  expect_true(is.numeric(res$p.value))
  expect_gte(res$p.value, 0)
  expect_lte(res$p.value, 1)
})
# -------------------------------------------------------------------------
# Known reference values
# -------------------------------------------------------------------------
test_that("tourists example matches Gibbons & Chakraborti (normal)", {
  
  tourists <- c(12362,12739,13057,13955,14123,15698,17523,18610,19842,
                20310,22500,23080,21916)
  
  res <- bartelsRankTest(tourists, alternative = "trend", method = "normal")
  
  expect_equal(unname(res$statistic["z"]), -3.6453, tolerance = 1e-3)
  expect_equal(res$p.value, 0.0001335203, tolerance = 1e-6)
  
})

test_that("tourists example matches Gibbons & Chakraborti (beta)", {
  
  tourists <- c(12362,12739,13057,13955,14123,15698,17523,18610,19842,
                20310,22500,23080,21916)
  
  res <- bartelsRankTest(tourists, alternative = "trend", method = "beta")
  
  expect_equal(unname(res$statistic["z"]), -3.6453, tolerance = 1e-3)
  expect_equal(res$p.value, 0.00000001209752, tolerance = 1e-6)
})

test_that("stock example gives correct RVN", {
  
  x <- c(528, 348, 264, -20, -167, 575, 410, -4, 430, -122)
  
  res <- bartelsRankTest(x, method = "beta")
  
  expect_equal(res$rvn, 2.048485, tolerance = 1e-4)
  expect_equal(unname(res$statistic["z"]), 0.083357, tolerance = 1e-4)
})
# -------------------------------------------------------------------------
# Methods agree on z-statistic
# -------------------------------------------------------------------------
test_that("normal and beta methods give same z-statistic", {
  
  set.seed(1)
  x <- rnorm(30)
  
  r_normal <- bartelsRankTest(x, method = "normal")
  r_beta   <- bartelsRankTest(x, method = "beta")
  
  expect_equal(
    unname(r_normal$statistic["z"]),
    unname(r_beta$statistic["z"]),
    tolerance = 1e-10
  )
})

test_that("auto selects beta for n <= 100", {
  
  set.seed(1)
  x <- rnorm(50)
  
  r_auto <- bartelsRankTest(x, method = "auto")
  r_beta <- bartelsRankTest(x, method = "beta")
  
  expect_equal(r_auto$p.value, r_beta$p.value)
})

test_that("auto selects normal for n > 100", {
  
  set.seed(1)
  x <- rnorm(150)
  
  r_auto   <- bartelsRankTest(x, method = "auto")
  r_normal <- bartelsRankTest(x, method = "normal")
  
  expect_equal(r_auto$p.value, r_normal$p.value)
})
# -------------------------------------------------------------------------
# Alternatives
# -------------------------------------------------------------------------
test_that("trend alternative detects increasing trend", {
  
  x <- c(1,2,3,4,5,6,7,8,9,10) + rnorm(10, sd = 0.1)
  
  res <- bartelsRankTest(x, alternative = "trend", method = "normal")
  
  expect_lt(res$p.value, 0.05)
})

test_that("oscillation alternative detects oscillation", {
  
  x <- rep(c(1, 10), 10) + rnorm(20, sd = 0.1)
  
  res <- bartelsRankTest(x, alternative = "oscillation", method = "normal")
  
  expect_lt(res$p.value, 0.05)
})

test_that("two-sided p-value >= one-sided p-value", {
  
  tourists <- c(12362,12739,13057,13955,14123,15698,17523,18610,19842,
                20310,22500,23080,21916)
  
  r_trend <- bartelsRankTest(tourists, alternative = "trend",     method = "normal")
  r_two   <- bartelsRankTest(tourists, alternative = "two.sided", method = "normal")
  
  expect_gte(r_two$p.value, r_trend$p.value)
})

test_that("trend and oscillation p-values sum to 1", {
  
  set.seed(1)
  x <- rnorm(20)
  
  r_trend <- bartelsRankTest(x, alternative = "trend",       method = "normal")
  r_osc   <- bartelsRankTest(x, alternative = "oscillation", method = "normal")
  
  expect_equal(r_trend$p.value + r_osc$p.value, 1, tolerance = 1e-10)
})

test_that("invalid alternative throws error", {
  
  expect_error(
    bartelsRankTest(rnorm(20), alternative = "invalid")
  )
})

test_that("invalid method throws error", {
  
  expect_error(
    bartelsRankTest(rnorm(20), method = "invalid")
  )
})
# -------------------------------------------------------------------------
# NA handling
# -------------------------------------------------------------------------
test_that("NA values are silently removed", {
  
  set.seed(1)
  x    <- rnorm(25)
  x_na <- x
  x_na[c(3, 15)] <- NA
  
  res_clean <- bartelsRankTest(x[!is.na(x_na)], method = "normal")
  res_na    <- bartelsRankTest(x_na,             method = "normal")
  
  expect_equal(unname(res_clean$statistic), unname(res_na$statistic))
  expect_equal(res_clean$p.value, res_na$p.value)
})
# -------------------------------------------------------------------------
# Input validation
# -------------------------------------------------------------------------
test_that("n < 10 throws error", {
  
  expect_error(
    bartelsRankTest(1:9),
    "greater than 9"
  )
})

test_that("parameter n is correct", {
  
  tourists <- c(12362,12739,13057,13955,14123,15698,17523,18610,19842,
                20310,22500,23080,21916)
  
  res <- bartelsRankTest(tourists)
  
  expect_equal(unname(res$parameter["n"]), 13L)
})

test_that("extra components present in output", {
  
  set.seed(1)
  x <- rnorm(20)
  
  res <- bartelsRankTest(x)
  
  expect_true(!is.null(res$rvn))
  expect_true(!is.null(res$nm))
  expect_true(!is.null(res$mu))
  expect_true(!is.null(res$var))
})
# -------------------------------------------------------------------------
# Matches randtests reference
# -------------------------------------------------------------------------
test_that("matches randtests::bartels.rank.test on z and p (normal)", {
  
  skip_if_not_installed("randtests")
  
  set.seed(1)
  x <- rnorm(30)
  
  res_ours <- bartelsRankTest(x, alternative = "two.sided", method = "normal")
  res_ref  <- randtests::bartels.rank.test(x, alternative = "two.sided",
                                           pvalue = "normal")
  
  expect_equal(
    unname(res_ours$statistic["z"]),
    unname(res_ref$statistic),
    tolerance = 1e-4
  )
  expect_equal(res_ours$p.value, res_ref$p.value, tolerance = 1e-6)
})
# -------------------------------------------------------------------------
# Print compatibility
# -------------------------------------------------------------------------
test_that("print.htest works", {
  
  set.seed(1)
  x <- rnorm(20)
  
  res <- bartelsRankTest(x)
  
  expect_output(print(res), "Bartels")
})