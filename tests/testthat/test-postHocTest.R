library(testthat)
library(lumen)

set.seed(1)
df <- data.frame(
  y = c(rnorm(20, 5), rnorm(20, 7), rnorm(20, 9)),
  g = factor(rep(c("A","B","C"), each = 20))
)
fit <- aov(y ~ g, data = df)

test_that("postHocTest: returns PostHocTest", {
  expect_s3_class(postHocTest(fit), "PostHocTest")
})

test_that("postHocTest: result is a list", {
  expect_true(is.list(postHocTest(fit)))
})

test_that("postHocTest: hsd method works", {
  res <- postHocTest(fit, method = "hsd")
  expect_s3_class(res, "PostHocTest")
})

test_that("postHocTest: bonferroni method works", {
  res <- postHocTest(fit, method = "bonferroni")
  expect_s3_class(res, "PostHocTest")
})

test_that("postHocTest: lsd method works", {
  res <- postHocTest(fit, method = "lsd")
  expect_s3_class(res, "PostHocTest")
})

test_that("postHocTest: scheffe method works", {
  res <- postHocTest(fit, method = "scheffe")
  expect_s3_class(res, "PostHocTest")
})

test_that("postHocTest: p-values in [0,1]", {
  res <- postHocTest(fit, method = "hsd")
  pvals <- res[[1]][, "pval"]
  expect_true(all(pvals >= 0 & pvals <= 1, na.rm = TRUE))
})

test_that("postHocTest: A vs C significant (large difference)", {
  res <- postHocTest(fit, method = "hsd")
  # A-C comparison should have small p-value given mean diff of ~4
  pvals <- res[[1]][, "pval"]
  ac_p  <- pvals[grep("A-C|C-A", names(pvals))]
  expect_lt(min(ac_p), 0.05)
})
