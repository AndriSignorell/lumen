library(testthat)
library(lumen)

# Agresti (2007) p.39 - gender x party
M <- as.table(rbind(c(762, 327, 468), c(484, 239, 477)))
dimnames(M) <- list(gender = c("M","F"), party = c("Democrat","Independent","Republican"))

test_that("gTest: returns htest", {
  expect_s3_class(gTest(M), "htest")
})

test_that("gTest: statistic named G", {
  expect_named(gTest(M)$statistic, "G")
})

test_that("gTest: p.value in [0,1]", {
  res <- gTest(M)
  expect_true(res$p.value >= 0 && res$p.value <= 1)
})

test_that("gTest: independence test df = (r-1)(c-1)", {
  res <- gTest(M)
  expect_equal(unname(res$parameter), (2-1) * (3-1))
})

test_that("gTest: Agresti example significant", {
  expect_lt(gTest(M)$p.value, 0.001)
})

test_that("gTest: GOF uniform p agrees with chisq.test for large n", {
  x <- c(100, 100, 100, 100)
  expect_gt(gTest(x)$p.value, 0.99)
})

test_that("gTest: GOF non-uniform gives small p", {
  x <- c(90, 10, 50, 50)
  expect_lt(gTest(x)$p.value, 0.05)
})

test_that("gTest: GOF df = k-1", {
  x <- c(20, 15, 25)
  expect_equal(unname(gTest(x)$parameter), 2L)
})

test_that("gTest: G statistic >= 0", {
  expect_gte(unname(gTest(M)$statistic), 0)
  expect_gte(unname(gTest(c(20, 15, 25))$statistic), 0)
})

test_that("gTest: observed and expected in result", {
  res <- gTest(M)
  expect_false(is.null(res$observed))
  expect_false(is.null(res$expected))
  expect_equal(sum(res$expected), sum(res$observed), tolerance = 1e-10)
})

test_that("gTest: rescaleP=TRUE allows non-summing p", {
  x <- c(89, 37, 30, 28, 2)
  p <- c(40, 20, 20, 15, 5)
  res <- gTest(x, p = p, rescaleP = TRUE)
  expect_s3_class(res, "htest")
})

test_that("gTest: williams correction gives larger p than none", {
  res_none <- gTest(M, correct = "none")
  res_will <- gTest(M, correct = "williams")
  expect_gte(res_will$p.value, res_none$p.value)
})
