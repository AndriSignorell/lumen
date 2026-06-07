library(testthat)
library(lumen)

test_that("pearsonTest: returns htest", {
  set.seed(1)
  expect_s3_class(pearsonTest(rnorm(100)), "htest")
})

test_that("pearsonTest: p.value in [0,1]", {
  set.seed(1)
  res <- pearsonTest(rnorm(100))
  expect_true(res$p.value >= 0 && res$p.value <= 1)
})

test_that("pearsonTest: normal data gives large p", {
  set.seed(42)
  expect_gt(pearsonTest(rnorm(300))$p.value, 0.05)
})

test_that("pearsonTest: non-normal data gives small p", {
  set.seed(1)
  expect_lt(pearsonTest(runif(300))$p.value, 0.05)
})

test_that("pearsonTest: statistic >= 0", {
  set.seed(1)
  expect_gte(unname(pearsonTest(rnorm(50))$statistic), 0)
})

test_that("pearsonTest: NAs handled via complete.cases", {
  set.seed(1)
  x <- c(rnorm(100), NA)
  expect_s3_class(pearsonTest(x), "htest")
})

test_that("pearsonTest: adjust=FALSE changes df", {
  set.seed(1)
  x  <- rnorm(100)
  r1 <- pearsonTest(x, adjust = TRUE)
  r2 <- pearsonTest(x, adjust = FALSE)
  # different df -> different p-values (usually)
  # adjust changes the effective df used internally, reflected in p.value not parameter
  expect_false(identical(r1$p.value, r2$p.value))
})
