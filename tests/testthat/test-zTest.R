library(testthat)
library(lumen)

tol <- 1e-4
set.seed(1); x <- rnorm(50, mean = 5, sd = 2)

test_that("zTest: one-sample returns htest", {
  expect_s3_class(zTest(x, sd_pop = 2), "htest")
})

test_that("zTest: p.value in [0,1]", {
  res <- zTest(x, sd_pop = 2)
  expect_true(res$p.value >= 0 && res$p.value <= 1)
})

test_that("zTest: H0 true (mu=5) gives large p", {
  set.seed(42); xh <- rnorm(200, mean = 5, sd = 2)
  expect_gt(zTest(xh, mu = 5, sd_pop = 2)$p.value, 0.05)
})

test_that("zTest: H0 false gives small p", {
  set.seed(42); xh <- rnorm(200, mean = 5, sd = 2)
  expect_lt(zTest(xh, mu = 0, sd_pop = 2)$p.value, 0.001)
})

test_that("zTest: statistic matches manual z calculation", {
  res <- zTest(x, mu = 0, sd_pop = 2)
  z_manual <- (mean(x) - 0) / (2 / sqrt(length(x)))
  expect_equal(unname(res$statistic), z_manual, tolerance = tol)
})

test_that("zTest: CI contains true mean for large n", {
  set.seed(1); xh <- rnorm(10000, mean = 5, sd = 2)
  res <- zTest(xh, mu = 5, sd_pop = 2)
  expect_true(res$conf.int[1] <= 5 && 5 <= res$conf.int[2])
})

test_that("zTest: alternative='less' p <= two.sided p", {
  # use mu > mean(x) so zstat < 0 and the left-tail p-value is small
  res2 <- zTest(x, mu = 10, sd_pop = 2, alternative = "two.sided")
  resL <- zTest(x, mu = 10, sd_pop = 2, alternative = "less")
  expect_lte(resL$p.value, res2$p.value + tol)
})

test_that("zTest: two-sample works", {
  set.seed(2); y <- rnorm(40, mean = 3, sd = 1.5)
  res <- zTest(x, y, sd_pop = 2)
  expect_s3_class(res, "htest")
})
