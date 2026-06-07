library(testthat)
library(lumen)

tol <- 1e-4
set.seed(1); x <- rnorm(30, sd = 2); y <- rnorm(25, sd = 3)

test_that("varTest: one-sample returns htest", {
  expect_s3_class(varTest(x, sigma2_0 = 4), "htest")
})

test_that("varTest: p.value in [0,1]", {
  res <- varTest(x, sigma2_0 = 4)
  expect_true(res$p.value >= 0 && res$p.value <= 1)
})

test_that("varTest: one-sample H0 true -> large p", {
  set.seed(42); xh <- rnorm(200, sd = 2)
  res <- varTest(xh, sigma2_0 = 4)
  expect_gt(res$p.value, 0.05)
})

test_that("varTest: one-sample H0 false -> small p", {
  set.seed(42); xh <- rnorm(200, sd = 5)
  res <- varTest(xh, sigma2_0 = 1)
  expect_lt(res$p.value, 0.05)
})

test_that("varTest: two-sample returns htest", {
  expect_s3_class(varTest(x, y), "htest")
})

test_that("varTest: two-sample equal variances gives large p", {
  set.seed(42); a <- rnorm(100, sd=2); b <- rnorm(100, sd=2)
  expect_gt(varTest(a, b)$p.value, 0.05)
})

test_that("varTest: two-sample unequal variances gives small p", {
  set.seed(42); a <- rnorm(200, sd=1); b <- rnorm(200, sd=5)
  expect_lt(varTest(a, b)$p.value, 0.05)
})

test_that("varTest: formula method works", {
  df <- data.frame(x = c(x, y), g = factor(rep(c("A","B"), c(30, 25))))
  expect_s3_class(varTest(x ~ g, data = df), "htest")
})

test_that("varTest: alternative='less' and 'greater' give valid p", {
  for (alt in c("less","greater")) {
    res <- varTest(x, sigma2_0 = 4, alternative = alt)
    expect_true(res$p.value >= 0 && res$p.value <= 1)
  }
})
