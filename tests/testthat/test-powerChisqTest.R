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

# powerChisqTest() -----------------------------------------------------------

pw <- function(w, n, df, a)
  pchisq(qchisq(a, df, lower.tail = FALSE), df, ncp = n * w^2, lower.tail = FALSE)

test_that("power by the noncentral chi-square", {
  r <- powerChisqTest(effectSize = 0.289, df = 3, n = 100, sig.level = 0.05)
  expect_s3_class(r, "power.htest")
  expect_equal(r$power, pw(0.289, 100, 3, 0.05))
  expect_identical(r$method, "Chi squared power calculation")
  expect_named(r, c("effectSize", "n", "df", "sig.level", "power",
                    "method", "note"))
  # w = 0 gives power = alpha
  expect_equal(powerChisqTest(effectSize = 0, n = 50, df = 2)$power, 0.05)
})

test_that("solving for n, effectSize and sig.level inverts the power", {
  n <- powerChisqTest(effectSize = 0.1, df = 20, power = 0.8)$n
  expect_equal(pw(0.1, n, 20, 0.05), 0.8, tolerance = 1e-4)
  
  w <- powerChisqTest(n = 140, df = 2, sig.level = 0.01, power = 0.9)$effectSize
  expect_equal(pw(w, 140, 2, 0.01), 0.9, tolerance = 1e-4)
  
  a <- powerChisqTest(effectSize = 0.3, n = 120, df = 2, power = 0.8,
                      sig.level = NULL)$sig.level
  expect_equal(pw(0.3, 120, 2, a), 0.8, tolerance = 1e-4)
})

test_that("argument checks", {
  expect_error(powerChisqTest(n = 10, effectSize = 0.3), "'df' must always")
  expect_error(powerChisqTest(df = 2, effectSize = 0.3), "exactly one")
  expect_error(powerChisqTest(n = 10, effectSize = 0.3, df = 2, power = 0.8),
               "exactly one")
  expect_error(powerChisqTest(n = 10, effectSize = -0.1, df = 2), "positive")
  expect_error(powerChisqTest(n = 0.5, effectSize = 0.1, df = 2), "at least 1")
  expect_error(powerChisqTest(n = 10, effectSize = 0.1, df = 2, sig.level = 1.5),
               "sig.level")
  expect_error(powerChisqTest(n = 10, effectSize = 0.1, df = 2, sig.level = "a"),
               "sig.level")
  expect_error(powerChisqTest(effectSize = 0.1, df = 2, power = 2), "power")
})
