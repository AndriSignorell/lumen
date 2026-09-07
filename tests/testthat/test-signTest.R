library(testthat)
library(lumen)

# ===============================================================
# signTest TESTS
# ===============================================================

# helper: the structural invariants every result has to satisfy.
# Expectations raised inside a helper are counted by testthat, so this
# stays a function rather than being inlined everywhere.
.check_signtest <- function(res) {

  expect_s3_class(res, "htest")
  expect_false(is.null(res$statistic))
  expect_false(is.null(res$p.value))
  expect_false(is.null(res$conf.int))
  expect_false(is.null(res$estimate))
  expect_true(res$p.value >= 0 && res$p.value <= 1)
  expect_true(res$conf.int[1] <= res$conf.int[2])

  invisible(TRUE)
}

# reference data from the documentation examples
x <- c(1.83,  0.50,  1.62,  2.48, 1.68, 1.88, 1.55, 3.06, 1.30)
y <- c(0.878, 0.647, 0.598, 2.05, 1.06, 1.29, 1.06, 3.14, 1.29)

d.light <- data.frame(
  black = c(25.85, 28.84, 32.05, 25.74, 20.89, 41.05, 25.01, 24.96, 27.47),
  white = c(18.23, 20.84, 22.96, 19.68, 19.5,  24.98, 16.61, 16.07, 24.59),
  d     = c(7.62, 8.00,  9.09,  6.06,  1.39,  16.07,  8.40,  8.89,  2.88)
)


test_that("signTest: two-sample result has the documented structure", {

  res <- signTest(x, y)

  .check_signtest(res)
  expect_equal(res$method, "Dependent-samples Sign-Test")
  expect_equal(names(res$statistic), "S")
  expect_equal(names(res$parameter), "number of differences")
})


test_that("signTest: one-sample result has the documented structure", {

  res <- signTest(x = d.light$d, mu = 4)

  .check_signtest(res)
  expect_equal(res$method, "One-sample Sign-Test")
})


test_that("signTest: the paired test equals the one-sample test on the differences", {

  res.diff <- signTest(x = d.light$black - d.light$white)
  res.pair <- signTest(x = d.light$black, y = d.light$white)

  expect_equal(res.diff$p.value, res.pair$p.value, tolerance = 1e-10)
  expect_equal(unname(res.diff$statistic), unname(res.pair$statistic),
               tolerance = 1e-10)
})


test_that("signTest: S counts the positive differences", {

  d <- d.light$d - 4
  expect_equal(unname(signTest(x = d.light$d, mu = 4)$statistic), sum(d > 0))

  # 3 positive, 2 negative
  expect_equal(unname(signTest(x = c(-2, -1, 1, 2, 3))$statistic), 3)

  # all positive, all negative
  expect_equal(unname(signTest(x = c(1, 2, 3, 4, 5))$statistic), 5)
  expect_equal(unname(signTest(x = c(-1, -2, -3, -4, -5))$statistic), 0)
})


test_that("signTest: the p-value matches binom.test on the sign counts", {

  d       <- x - y
  n.valid <- sum(d != 0)
  s       <- sum(d > 0)
  ref     <- binom.test(x = s, n = n.valid, p = 0.5)

  expect_equal(signTest(x, y)$p.value, ref$p.value, tolerance = 1e-10)
})


test_that("signTest: the one-sided alternatives are consistent", {

  res.ts <- signTest(x = d.light$d, mu = 4, alternative = "two.sided")
  res.gt <- signTest(x = d.light$d, mu = 4, alternative = "greater")
  res.lt <- signTest(x = d.light$d, mu = 4, alternative = "less")

  .check_signtest(res.ts)
  .check_signtest(res.gt)
  .check_signtest(res.lt)

  # the two one-sided p-values overlap at the observed value, so they sum
  # to more than one (the standard relation for a discrete statistic)
  expect_gt(res.gt$p.value + res.lt$p.value, 1)

  # at least one of them is no larger than the two-sided p-value
  expect_true(res.gt$p.value <= res.ts$p.value ||
                res.lt$p.value <= res.ts$p.value)
})


test_that("signTest: a higher conf.level gives a wider interval", {

  res.95 <- signTest(x = d.light$d, mu = 4, conf.level = 0.95)
  res.99 <- signTest(x = d.light$d, mu = 4, conf.level = 0.99)

  expect_lte(res.99$conf.int[1], res.95$conf.int[1])
  expect_gte(res.99$conf.int[2], res.95$conf.int[2])
})


test_that("signTest: values tied with mu drop out of n", {

  # x = c(1, 2, 4) with mu = 2: d = c(-1, 0, 2), two usable differences
  res <- signTest(x = c(1, 2, 4), mu = 2)

  expect_equal(unname(res$parameter), 2)
  expect_equal(unname(res$statistic), 1)
})


test_that("signTest: missing values are removed", {

  res.na    <- signTest(x = c(1, 2, NA, 3, 4))
  res.clean <- signTest(x = c(1, 2, 3, 4))

  expect_equal(res.na$p.value, res.clean$p.value, tolerance = 1e-10)
})


test_that("signTest: invalid input throws", {

  expect_error(signTest(x = c(1, 2, 3), mu = c(1, 2)))
  expect_error(signTest(x = c(1, 2, 3), conf.level = 1.5))
  expect_error(signTest(x = c(1, 2, 3), y = c(1, 2)))
  expect_error(signTest(x = "a"))
})


test_that("signTest: null.value carries mu and is labelled by design", {

  res <- signTest(x = d.light$d, mu = 4)
  expect_equal(unname(res$null.value), 4)
  expect_equal(names(res$null.value), "median")

  expect_equal(names(signTest(x, y)$null.value), "median difference")
})


test_that("signTest: the estimate is the median of x for mu = 0", {

  res <- signTest(x = c(1, 2, 3, 4, 5))
  expect_equal(unname(res$estimate), median(c(1, 2, 3, 4, 5)),
               tolerance = 1e-10)
})
