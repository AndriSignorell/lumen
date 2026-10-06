library(testthat)
library(lumen)

test_that("breuschPaganTest: returns htest", {
  fit <- lm(Sepal.Length ~ Sepal.Width, data = iris)
  expect_s3_class(breuschPaganTest(fit), "htest")
})

test_that("breuschPaganTest: result has statistic, parameter, p.value", {
  fit <- lm(Sepal.Length ~ Sepal.Width, data = iris)
  res <- breuschPaganTest(fit)
  expect_named(res$statistic, "BP")
  expect_named(res$parameter, "df")
  expect_false(is.null(res$p.value))
})

test_that("breuschPaganTest: p.value in [0,1]", {
  fit <- lm(Sepal.Length ~ Sepal.Width, data = iris)
  res <- breuschPaganTest(fit)
  expect_true(res$p.value >= 0 && res$p.value <= 1)
})

test_that("breuschPaganTest: homoscedastic data gives large p-value", {
  set.seed(1)
  x  <- 1:100
  y  <- 2 + 3 * x + rnorm(100, sd = 1)
  res <- breuschPaganTest(lm(y ~ x))
  expect_gt(res$p.value, 0.05)
})

test_that("breuschPaganTest: heteroscedastic data gives small p-value", {
  set.seed(1)
  x  <- 1:200
  y  <- 1 + x + rnorm(200, sd = 0.1 * x)
  res <- breuschPaganTest(lm(y ~ x))
  expect_lt(res$p.value, 0.05)
})

test_that("breuschPaganTest: df = number of predictors", {
  fit <- lm(Sepal.Length ~ Sepal.Width + Petal.Length, data = iris)
  res <- breuschPaganTest(fit)
  # auxiliary regression: resid^2 ~ fitted -> 1 predictor -> df=1
  expect_equal(unname(res$parameter["df"]), 2L)
})

test_that("breuschPaganTest: non-lm input throws error", {
  expect_error(breuschPaganTest(list(a = 1)))
})

test_that("breuschPaganTest: BP statistic is non-negative", {
  fit <- lm(mpg ~ wt + hp, data = mtcars)
  expect_gte(unname(breuschPaganTest(fit)$statistic), 0)
})
