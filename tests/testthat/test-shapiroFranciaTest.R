library(testthat)
library(lumen)

test_that("shapiroFranciaTest: returns htest", {
  set.seed(1); expect_s3_class(shapiroFranciaTest(rnorm(50)), "htest")
})

test_that("shapiroFranciaTest: p.value in [0,1]", {
  set.seed(1); res <- shapiroFranciaTest(rnorm(50))
  expect_true(res$p.value >= 0 && res$p.value <= 1)
})

test_that("shapiroFranciaTest: normal data gives large p", {
  set.seed(42); expect_gt(shapiroFranciaTest(rnorm(200))$p.value, 0.05)
})

test_that("shapiroFranciaTest: non-normal data gives small p", {
  set.seed(1); expect_lt(shapiroFranciaTest(runif(200))$p.value, 0.05)
})

test_that("shapiroFranciaTest: statistic W in (0,1]", {
  set.seed(1); res <- shapiroFranciaTest(rnorm(50))
  expect_gt(unname(res$statistic), 0)
  expect_lte(unname(res$statistic), 1)
})

test_that("shapiroFranciaTest: n < 5 throws error", {
  expect_error(shapiroFranciaTest(rnorm(4)))
})

test_that("shapiroFranciaTest: n > 5000 throws error", {
  expect_error(shapiroFranciaTest(rnorm(5001)))
})

test_that("shapiroFranciaTest: NAs silently removed", {
  set.seed(1); x <- c(rnorm(50), NA)
  expect_s3_class(shapiroFranciaTest(x), "htest")
})
