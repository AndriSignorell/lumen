library(testthat)
library(lumen)

set.seed(42)
x    <- rep(c(-1, 1), 50)
err1 <- rnorm(100)
err2 <- stats::filter(err1, 0.9, method = "recursive")
df1  <- data.frame(y = 1 + x + err1, x = x)
df2  <- data.frame(y = 1 + x + err2, x = x)

test_that("durbinWatsonTest: returns htest (formula method)", {
  res <- durbinWatsonTest(y ~ x, data = df1)
  expect_s3_class(res, "htest")
})

test_that("durbinWatsonTest: DW statistic in (0, 4)", {
  res <- durbinWatsonTest(y ~ x, data = df1)
  expect_gt(unname(res$statistic), 0)
  expect_lt(unname(res$statistic), 4)
})

test_that("durbinWatsonTest: p.value in [0,1]", {
  res <- durbinWatsonTest(y ~ x, data = df1)
  expect_true(res$p.value >= 0 && res$p.value <= 1)
})

test_that("durbinWatsonTest: iid errors -> large p-value (greater)", {
  res <- durbinWatsonTest(y ~ x, data = df1, alternative = "greater")
  expect_gt(res$p.value, 0.05)
})

test_that("durbinWatsonTest: AR(1) errors -> small p-value (greater)", {
  res <- durbinWatsonTest(y ~ x, data = df2, alternative = "greater")
  expect_lt(res$p.value, 0.05)
})

test_that("durbinWatsonTest: DW ~ 2 for iid errors", {
  # DW ≈ 2*(1-rho), iid -> rho≈0 -> DW≈2
  res <- durbinWatsonTest(y ~ x, data = df1)
  expect_equal(unname(res$statistic), 2, tolerance = 0.5)
})

test_that("durbinWatsonTest: DW < 2 for positive autocorrelation", {
  res <- durbinWatsonTest(y ~ x, data = df2)
  expect_lt(unname(res$statistic), 2)
})

test_that("durbinWatsonTest: lm method agrees with formula method", {
  fit <- lm(y ~ x, data = df1)
  res_form <- durbinWatsonTest(y ~ x, data = df1)
  res_lm   <- durbinWatsonTest(fit)
  expect_equal(unname(res_form$statistic), unname(res_lm$statistic), tolerance = 1e-10)
})

test_that("durbinWatsonTest: numeric method works", {
  e <- residuals(lm(y ~ x, data = df1))
  res <- durbinWatsonTest(e)
  expect_s3_class(res, "htest")
})

test_that("durbinWatsonTest: two.sided p-value <= 2 * one-sided", {
  res_two <- durbinWatsonTest(y ~ x, data = df1, alternative = "two.sided")
  res_one <- durbinWatsonTest(y ~ x, data = df1, alternative = "greater")
  expect_lte(res_two$p.value, 2 * res_one$p.value + 1e-10)
})

test_that("durbinWatsonTest: unsupported class throws error", {
  expect_error(durbinWatsonTest(list(a = 1)))
})
