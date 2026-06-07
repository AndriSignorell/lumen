library(testthat)
library(lumen)

# Fixed test data: iid errors -> no serial correlation
set.seed(42)
x  <- rep(c(1, -1), 50)
y1 <- 1 + x + rnorm(100)
# AR(1) errors: strong autocorrelation
y2 <- stats::filter(y1, 0.9, method = "recursive")

test_that("breuschGodfreyTest: returns htest", {
  res <- breuschGodfreyTest(y1 ~ x)
  expect_s3_class(res, "htest")
})

test_that("breuschGodfreyTest: result components present", {
  res <- breuschGodfreyTest(y1 ~ x)
  expect_false(is.null(res$statistic))
  expect_false(is.null(res$p.value))
  expect_false(is.null(res$parameter))
})

test_that("breuschGodfreyTest: p.value in [0,1]", {
  res <- breuschGodfreyTest(y1 ~ x)
  expect_true(res$p.value >= 0 && res$p.value <= 1)
})

test_that("breuschGodfreyTest: iid residuals give large p-value", {
  res <- breuschGodfreyTest(y1 ~ x, order = 1)
  expect_gt(res$p.value, 0.05)
})

test_that("breuschGodfreyTest: AR(1) residuals give small p-value", {
  res <- breuschGodfreyTest(y2 ~ x, order = 1)
  expect_lt(res$p.value, 0.05)
})

test_that("breuschGodfreyTest: order=4 gives df=4", {
  res <- breuschGodfreyTest(y1 ~ x, order = 4)
  expect_equal(unname(res$parameter), 4L)
})

test_that("breuschGodfreyTest: type='F' returns F statistic", {
  res <- breuschGodfreyTest(y1 ~ x, type = "F")
  expect_named(res$parameter, c("df1", "df2"))
})

test_that("breuschGodfreyTest: accepts lm object", {
  fit <- lm(y1 ~ x)
  res <- breuschGodfreyTest(fit)
  expect_s3_class(res, "htest")
})

test_that("breuschGodfreyTest: coefficients and vcov present", {
  res <- breuschGodfreyTest(y1 ~ x)
  expect_false(is.null(res$coefficients))
  expect_false(is.null(res$vcov))
})

test_that("breuschGodfreyTest: statistic >= 0", {
  res <- breuschGodfreyTest(y1 ~ x)
  expect_gte(unname(res$statistic), 0)
})
