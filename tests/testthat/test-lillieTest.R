library(testthat)
library(lumen)

test_that("lillieTest: returns htest", {
  expect_s3_class(lillieTest(rnorm(50)), "htest")
})

test_that("lillieTest: p.value in [0,1]", {
  set.seed(1)
  res <- lillieTest(rnorm(50))
  expect_true(res$p.value >= 0 && res$p.value <= 1)
})

test_that("lillieTest: normal data gives large p", {
  set.seed(42)
  expect_gt(lillieTest(rnorm(200))$p.value, 0.05)
})

test_that("lillieTest: non-normal data gives small p", {
  set.seed(1)
  expect_lt(lillieTest(runif(200))$p.value, 0.05)
})

test_that("lillieTest: statistic named K", {
  res <- lillieTest(rnorm(20))
  expect_named(res$statistic, "D")
})

test_that("lillieTest: statistic > 0", {
  set.seed(1)
  expect_gt(lillieTest(rnorm(30))$statistic, 0)
})

test_that("lillieTest: n < 5 throws error", {
  expect_error(lillieTest(rnorm(4)))
})

test_that("lillieTest: NAs silently removed", {
  set.seed(1)
  x <- c(rnorm(50), NA, NA)
  expect_s3_class(lillieTest(x), "htest")
})
