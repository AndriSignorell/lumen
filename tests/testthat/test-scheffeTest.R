library(testthat)
library(lumen)

set.seed(1)
df <- data.frame(
  y = c(rnorm(20, 5), rnorm(20, 7), rnorm(20, 9)),
  g = factor(rep(c("A","B","C"), each = 20))
)
fit <- aov(y ~ g, data = df)

test_that("scheffeTest: aov method returns PostHocTest", {
  expect_s3_class(scheffeTest(fit), "PostHocTest")
})

test_that("scheffeTest: formula method works", {
  res <- scheffeTest(y ~ g, data = df)
  expect_s3_class(res, "PostHocTest")
})

test_that("scheffeTest: p-values in [0,1]", {
  res  <- scheffeTest(fit)
  pvals <- res[[1]][, "pval"]
  expect_true(all(pvals >= 0 & pvals <= 1, na.rm = TRUE))
})

test_that("scheffeTest: A vs C significant", {
  res   <- scheffeTest(fit)
  pvals <- res[[1]][, "pval"]
  ac_p  <- pvals[grep("A-C|C-A", names(pvals))]
  expect_lt(min(ac_p), 0.05)
})

test_that("scheffeTest: CI contains diff for equal groups", {
  set.seed(42)
  df2 <- data.frame(y = rnorm(60), g = factor(rep(1:3, each=20)))
  res <- scheffeTest(aov(y ~ g, data = df2))
  # CI for all pairs should contain 0
  diffs <- res[[1]][, "diff"]
  lcis  <- res[[1]][, "lwr.ci"]
  ucis  <- res[[1]][, "upr.ci"]
  expect_true(all(lcis <= 0 & 0 <= ucis, na.rm = TRUE))
})

test_that("scheffeTest: wider CI with higher conf.level", {
  res95 <- scheffeTest(fit, conf.level = 0.95)[[1]]
  res99 <- scheffeTest(fit, conf.level = 0.99)[[1]]
  w95   <- mean(res95[,"upr.ci"] - res95[,"lwr.ci"], na.rm = TRUE)
  w99   <- mean(res99[,"upr.ci"] - res99[,"lwr.ci"], na.rm = TRUE)
  expect_gt(w99, w95)
})
