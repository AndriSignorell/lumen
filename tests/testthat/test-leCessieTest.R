library(testthat)
library(lumen)

set.seed(111)
x1  <- factor(sample(1:3, 50, replace = TRUE))
x2  <- rnorm(50)
obs <- sample(c(0, 1), 50, replace = TRUE)
fit <- glm(obs ~ x1 + x2, family = binomial)
f   <- fitted(fit)
X   <- model.matrix(fit)[, -1, drop = FALSE]

test_that("leCessieTest: returns htest / LeCessieTest", {
  res <- leCessieTest(fit = f, obs = obs, X = X)
  expect_s3_class(res, "htest")
  expect_s3_class(res, "LeCessieTest")
})

test_that("leCessieTest: statistic named Z", {
  res <- leCessieTest(fit = f, obs = obs, X = X)
  expect_named(res$statistic, "Z")
})

test_that("leCessieTest: p.value in [0, 1]", {
  res <- leCessieTest(fit = f, obs = obs, X = X)
  expect_gte(res$p.value, 0)
  expect_lte(res$p.value, 1)
})

test_that("leCessieTest: sse, expected, sd are positive scalars", {
  res <- leCessieTest(fit = f, obs = obs, X = X)
  expect_gt(res$sse,      0)
  expect_gt(res$expected, 0)
  expect_gt(res$sd,       0)
})

test_that("leCessieTest: well-specified model gives large p", {
  set.seed(1)
  n   <- 500
  x   <- rnorm(n)
  eta <- -1 + 2 * x
  y   <- rbinom(n, 1, plogis(eta))
  g   <- glm(y ~ x, family = binomial)
  res <- leCessieTest(
    fit = fitted(g), obs = y,
    X   = model.matrix(g)[, -1, drop = FALSE]
  )
  expect_gt(res$p.value, 0.05)
})

test_that("leCessieTest: misspecified model gives small p", {
  set.seed(42)
  n   <- 1000
  x1  <- rnorm(n)
  x2  <- rnorm(n)
  # true model has interaction, fitted model misses it
  eta <- -1 + 2 * x1 + 2 * x2 + 5 * x1 * x2
  y   <- rbinom(n, 1, plogis(eta))
  g   <- glm(y ~ x1 + x2, family = binomial)   # missing interaction
  res <- leCessieTest(
    fit = fitted(g), obs = y,
    X   = model.matrix(g)[, -1, drop = FALSE]
  )
  expect_lt(res$p.value, 0.05)
})

test_that("leCessieTest: input validation - length mismatch", {
  expect_error(leCessieTest(f[-1], obs, X), "same length")
})

test_that("leCessieTest: input validation - fit out of [0,1]", {
  bad <- f; bad[1] <- 1.5
  expect_error(leCessieTest(bad, obs, X), "probabilities")
})

test_that("leCessieTest: input validation - non-binary obs", {
  bad <- obs; bad[1] <- 2
  expect_error(leCessieTest(f, bad, X), "binary")
})

test_that("leCessieTest: input validation - X row mismatch", {
  expect_error(leCessieTest(f, obs, X[-1, ]), "same number of rows")
})

test_that("leCessieTest: print method runs without error", {
  res <- leCessieTest(fit = f, obs = obs, X = X)
  expect_output(print(res))
})
