library(testthat)
library(lumen)

test_that("cramerVonMisesTest: returns htest", {
  set.seed(1)
  res <- cramerVonMisesTest(rnorm(50))
  expect_s3_class(res, "htest")
})

test_that("cramerVonMisesTest: result has statistic and p.value", {
  set.seed(1)
  res <- cramerVonMisesTest(rnorm(50))
  expect_true(!is.null(res$statistic))
  expect_true(!is.null(res$p.value))
})

test_that("cramerVonMisesTest: p.value in [0,1]", {
  set.seed(1)
  res <- cramerVonMisesTest(rnorm(50))
  expect_true(res$p.value >= 0 && res$p.value <= 1)
})

test_that("cramerVonMisesTest: normal data gives large p-value", {
  set.seed(42)
  res <- cramerVonMisesTest(rnorm(200))
  expect_gt(res$p.value, 0.05)
})

test_that("cramerVonMisesTest: uniform data gives small p-value", {
  set.seed(1)
  res <- cramerVonMisesTest(runif(200))
  expect_lt(res$p.value, 0.05)
})

test_that("cramerVonMisesTest: statistic > 0", {
  set.seed(1)
  expect_gt(cramerVonMisesTest(rnorm(30))$statistic, 0)
})

test_that("cramerVonMisesTest: n < 8 throws error", {
  expect_error(cramerVonMisesTest(rnorm(7)))
})

test_that("cramerVonMisesTest: NAs are removed", {
  set.seed(1)
  x <- c(rnorm(50), NA, NA)
  # should not error, NAs silently dropped
  res <- cramerVonMisesTest(x)
  expect_s3_class(res, "htest")
})

test_that("cramerVonMisesTest: method string correct", {
  res <- cramerVonMisesTest(rnorm(20))
  expect_equal(res$method, "Cramer-von Mises normality test")
})
