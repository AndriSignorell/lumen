library(testthat)
library(lumen)

test_that("powerChisqTest: returns power.htest", {
  res <- powerChisqTest(n = 100, effectSize = 0.3, df = 2)
  expect_s3_class(res, "power.htest")
})

test_that("powerChisqTest: power in (0,1]", {
  res <- powerChisqTest(n = 100, effectSize = 0.3, df = 2)
  expect_true(res$power > 0 && res$power <= 1)
})

test_that("powerChisqTest: larger n gives more power", {
  p50  <- powerChisqTest(n =  50, effectSize = 0.3, df = 2)$power
  p200 <- powerChisqTest(n = 200, effectSize = 0.3, df = 2)$power
  expect_gt(p200, p50)
})

test_that("powerChisqTest: larger effect (effectSize) gives more power", {
  p1 <- powerChisqTest(n = 100, effectSize = 0.1, df = 2)$power
  p3 <- powerChisqTest(n = 100, effectSize = 0.3, df = 2)$power
  expect_gt(p3, p1)
})

test_that("powerChisqTest: solves for n", {
  res <- powerChisqTest(effectSize = 0.3, df = 2, power = 0.80)
  expect_false(is.null(res$n))
  expect_gt(res$n, 0)
})

test_that("powerChisqTest: solves for effectSize", {
  res <- powerChisqTest(n = 100, df = 2, power = 0.80)
  expect_false(is.null(res$effectSize))
  expect_gt(res$effectSize, 0)
})

test_that("powerChisqTest: two NULLs throws error", {
  expect_error(powerChisqTest(df = 2, sig.level = 0.05))
})

test_that("powerChisqTest: sig.level increase gives more power", {
  p05 <- powerChisqTest(n = 100, effectSize = 0.3, df = 2, sig.level = 0.05)$power
  p10 <- powerChisqTest(n = 100, effectSize = 0.3, df = 2, sig.level = 0.10)$power
  expect_gt(p10, p05)
})
