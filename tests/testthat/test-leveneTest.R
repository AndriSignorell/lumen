library(testthat)
library(lumen)

# Equal variances: homoscedastic groups
set.seed(1)
g1 <- rnorm(30, 0, 1); g2 <- rnorm(30, 5, 1); g3 <- rnorm(30, 10, 1)
df_eq <- data.frame(x = c(g1, g2, g3), g = factor(rep(1:3, each = 30)))

# Unequal variances: heteroscedastic groups
g4 <- rnorm(30, 0, 1); g5 <- rnorm(30, 0, 5)
df_uneq <- data.frame(x = c(g4, g5), g = factor(rep(1:2, each = 30)))

test_that("leveneTest: returns htest (formula method)", {
  expect_s3_class(leveneTest(x ~ g, data = df_eq), "htest")
})

test_that("leveneTest: p.value in [0,1]", {
  res <- leveneTest(x ~ g, data = df_eq)
  expect_true(res$p.value >= 0 && res$p.value <= 1)
})

test_that("leveneTest: equal variances gives large p", {
  set.seed(42)
  g1 <- rnorm(50, sd=1); g2 <- rnorm(50, sd=1)
  df <- data.frame(x = c(g1,g2), g = factor(rep(1:2, each=50)))
  expect_gt(leveneTest(x ~ g, data = df)$p.value, 0.05)
})

test_that("leveneTest: unequal variances gives small p", {
  set.seed(42)
  g1 <- rnorm(100, sd=1); g2 <- rnorm(100, sd=10)
  df <- data.frame(x = c(g1,g2), g = factor(rep(1:2, each=100)))
  expect_lt(leveneTest(x ~ g, data = df)$p.value, 0.05)
})

test_that("leveneTest: center=mean gives original Levene test", {
  res <- leveneTest(x ~ g, data = df_eq, center = mean)
  expect_s3_class(res, "htest")
})

test_that("leveneTest: default method (formula) works", {
  res <- leveneTest(x ~ g, data = df_uneq)
  expect_false(is.null(res$statistic))
})

test_that("leveneTest: default method (default) works with vector + factor", {
  res <- leveneTest(df_eq$x, df_eq$g)
  expect_s3_class(res, "htest")
})
