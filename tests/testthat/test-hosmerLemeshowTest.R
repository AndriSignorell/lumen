library(testthat)
library(lumen)

set.seed(111)
x1  <- factor(sample(1:3, 50, replace = TRUE))
x2  <- rnorm(50)
obs <- sample(c(0, 1), 50, replace = TRUE)
fit <- glm(obs ~ x1 + x2, family = binomial)
f   <- fitted(fit)

test_that("hosmerLemeshowTest: returns htest / HosmerLemeshowTest (type C)", {
  res <- hosmerLemeshowTest(x = f, obs = obs, type = "C")
  expect_s3_class(res, "htest")
  expect_s3_class(res, "HosmerLemeshowTest")
})

test_that("hosmerLemeshowTest: returns htest / HosmerLemeshowTest (type H)", {
  res <- hosmerLemeshowTest(x = f, obs = obs, type = "H")
  expect_s3_class(res, "htest")
  expect_s3_class(res, "HosmerLemeshowTest")
})

test_that("hosmerLemeshowTest: statistic named X-squared", {
  res <- hosmerLemeshowTest(x = f, obs = obs)
  expect_named(res$statistic, "X-squared")
})

test_that("hosmerLemeshowTest: parameter named df equals nGroups - 2", {
  res <- hosmerLemeshowTest(x = f, obs = obs, nGroups = 10)
  expect_named(res$parameter, "df")
  expect_equal(unname(res$parameter), res$nGroups - 2L)
})

test_that("hosmerLemeshowTest: p.value in [0, 1]", {
  res <- hosmerLemeshowTest(x = f, obs = obs)
  expect_gte(res$p.value, 0)
  expect_lte(res$p.value, 1)
})

test_that("hosmerLemeshowTest: observed and expected are matrices", {
  res <- hosmerLemeshowTest(x = f, obs = obs)
  expect_true(is.matrix(res$observed))
  expect_true(is.matrix(res$expected))
  expect_equal(colnames(res$observed), c("0s", "1s"))
  expect_equal(colnames(res$expected), c("0s", "1s"))
})

test_that("hosmerLemeshowTest: well-specified model gives large p", {
  set.seed(1)
  n   <- 500
  x   <- rnorm(n)
  eta <- -1 + 2 * x
  y   <- rbinom(n, 1, plogis(eta))
  g   <- glm(y ~ x, family = binomial)
  res <- hosmerLemeshowTest(x = fitted(g), obs = y)
  expect_gt(res$p.value, 0.05)
})

test_that("hosmerLemeshowTest: nGroups respected", {
  res <- hosmerLemeshowTest(x = f, obs = obs, nGroups = 5)
  expect_lte(res$nGroups, 5L)
  expect_gte(res$nGroups, 3L)
})

test_that("hosmerLemeshowTest: input validation - length mismatch", {
  expect_error(hosmerLemeshowTest(f[-1], obs), "same length")
})

test_that("hosmerLemeshowTest: input validation - fit out of [0,1]", {
  bad <- f; bad[1] <- -0.1
  expect_error(hosmerLemeshowTest(bad, obs), "probabilities")
})

test_that("hosmerLemeshowTest: input validation - non-binary obs", {
  bad <- obs; bad[1] <- 3
  expect_error(hosmerLemeshowTest(f, bad), "binary")
})

test_that("hosmerLemeshowTest: input validation - nGroups < 3", {
  expect_error(hosmerLemeshowTest(f, obs, nGroups = 2), "nGroups")
})

test_that("hosmerLemeshowTest: print method runs without error", {
  res <- hosmerLemeshowTest(x = f, obs = obs)
  expect_output(print(res))
})

test_that("hosmerLemeshowTest: print with details runs without error", {
  res <- hosmerLemeshowTest(x = f, obs = obs)
  expect_output(print(res, details = TRUE))
})


test_that("hosmerLemeshowTest.glm: matches the default method on the same data", {
  res_default <- hosmerLemeshowTest(x = f, obs = obs, type = "C")
  res_glm     <- hosmerLemeshowTest(fit)

  expect_equal(unname(res_glm$statistic), unname(res_default$statistic))
  expect_equal(res_glm$p.value, res_default$p.value)
})

test_that("hosmerLemeshowTest.glm: data.name reflects the model formula", {
  res <- hosmerLemeshowTest(fit)
  expect_equal(res$data.name, "obs ~ x1 + x2")
})

test_that("hosmerLemeshowTest.glm: rejects a non-binomial glm", {
  g <- glm(x2 ~ x1, family = gaussian)
  expect_error(hosmerLemeshowTest(g), "binomial")
})

test_that("hosmerLemeshowTest: type = 'H' uses fixed [0,1] bins and can drop empty groups", {
  set.seed(4)
  fit_narrow <- runif(80, 0.35, 0.65)   # occupies only a few of 10 [0,1] bins
  obs_narrow <- rbinom(80, 1, fit_narrow)

  expect_warning(
    res <- hosmerLemeshowTest(fit_narrow, obs_narrow, type = "H"),
    "empty group"
  )
  expect_lt(res$nGroups, 10L)
  expect_equal(unname(res$parameter), res$nGroups - 2L)
})
