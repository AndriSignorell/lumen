library(testthat)
library(lumen)

tol <- 1e-4

# One-sample: matches t.test
set.seed(1); x <- rnorm(30, mean = 5, sd = 2)
ref1 <- t.test(x)

test_that("tTestA: one-sample returns htest", {
  expect_s3_class(tTestA(mx = mean(x), sx = sd(x), nx = length(x)), "htest")
})

test_that("tTestA: one-sample statistic matches t.test", {
  res <- tTestA(mx = mean(x), sx = sd(x), nx = length(x))
  expect_equal(unname(res$statistic), unname(ref1$statistic), tolerance = tol)
})

test_that("tTestA: one-sample p.value matches t.test", {
  res <- tTestA(mx = mean(x), sx = sd(x), nx = length(x))
  expect_equal(res$p.value, ref1$p.value, tolerance = tol)
})

test_that("tTestA: one-sample CI matches t.test", {
  res <- tTestA(mx = mean(x), sx = sd(x), nx = length(x))
  expect_equal(unname(res$conf.int[1]), ref1$conf.int[1], tolerance = tol)
  expect_equal(unname(res$conf.int[2]), ref1$conf.int[2], tolerance = tol)
})

# Two-sample
set.seed(2); y <- rnorm(25, mean = 3, sd = 1.5)
ref2 <- t.test(x, y)

test_that("tTestA: two-sample returns htest", {
  res <- tTestA(mean(x), sd(x), length(x), mean(y), sd(y), length(y))
  expect_s3_class(res, "htest")
})

test_that("tTestA: two-sample p.value matches t.test (Welch)", {
  res <- tTestA(mean(x), sd(x), length(x), mean(y), sd(y), length(y))
  expect_equal(res$p.value, ref2$p.value, tolerance = tol)
})

test_that("tTestA: var.equal=TRUE matches t.test", {
  ref_eq <- t.test(x, y, var.equal = TRUE)
  res    <- tTestA(mean(x), sd(x), length(x), mean(y), sd(y), length(y),
                   var.equal = TRUE)
  expect_equal(res$p.value, ref_eq$p.value, tolerance = tol)
})

test_that("tTestA: alternative='less' gives p < two.sided p", {
  # use mu > mean(x) so tstat < 0 and the left-tail p-value is small
  res_two  <- tTestA(mean(x), sd(x), length(x), mu = 10, alternative = "two.sided")
  res_less <- tTestA(mean(x), sd(x), length(x), mu = 10, alternative = "less")
  expect_lte(res_less$p.value, res_two$p.value)
})

test_that("tTestA: mu shifts test statistic", {
  res0 <- tTestA(mean(x), sd(x), length(x), mu = 0)
  res5 <- tTestA(mean(x), sd(x), length(x), mu = 5)
  expect_false(isTRUE(all.equal(res0$statistic, res5$statistic)))
})

test_that("tTestA: paired not supported -> error", {
  expect_error(tTestA(mean(x), sd(x), length(x), paired = TRUE))
})
