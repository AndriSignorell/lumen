library(testthat)
library(lumen)

tol <- 1e-6

test_that("drevweibull: density >= 0", {
  expect_true(all(drevweibull(seq(-5, -0.01, by = 0.1), 0, 1, 1) >= 0))
})

test_that("drevweibull: density = 0 for x >= loc", {
  expect_equal(drevweibull(c(0, 1, 2), loc = 0), c(0, 0, 0))
})

test_that("drevweibull: integrates to 1", {
  x <- seq(-30, 0, length.out = 100001)
  dx <- x[2] - x[1]
  expect_equal(sum(drevweibull(x, 0, 1, 1)) * dx, 1, tolerance = 1e-3)
})

test_that("drevweibull log=TRUE", {
  x <- c(-3, -2, -1)
  expect_equal(drevweibull(x, -1, 0.5, 0.8, log = TRUE),
               log(drevweibull(x, -1, 0.5, 0.8)), tolerance = tol)
})

test_that("drevweibull: invalid scale/shape throws error", {
  expect_error(drevweibull(-1, scale = 0))
  expect_error(drevweibull(-1, shape = 0))
})

test_that("prevweibull: in [0,1] and non-decreasing", {
  q <- seq(-10, 0, by = 0.5)
  p <- prevweibull(q, 0, 1, 1)
  expect_true(all(p >= 0 & p <= 1))
  expect_true(all(diff(p) >= 0))
})

test_that("prevweibull: lower.tail=FALSE complement", {
  q <- c(-3, -2, -1)
  expect_equal(prevweibull(q, -1, 0.5, 0.8, lower.tail = FALSE),
               1 - prevweibull(q, -1, 0.5, 0.8), tolerance = tol)
})

test_that("prevweibull(qrevweibull(p)) roundtrip", {
  p <- c(0.1, 0.25, 0.5, 0.75, 0.9)
  expect_equal(prevweibull(qrevweibull(p, -1, 2, 0.8), -1, 2, 0.8), p, tolerance = tol)
})

test_that("qrevweibull: p = 0 and p = 1 give the end points of the support", {
  # support is (-Inf, loc)
  expect_equal(qrevweibull(c(0, 1), loc = -1), c(-Inf, -1))
})

test_that("qrevweibull: p outside [0,1] gives NaN with a warning", {
  expect_warning(res <- qrevweibull(-0.1), "NaN")
  expect_true(is.nan(res))
})

test_that("qrevweibull: log.p and lower.tail", {
  p <- c(0.1, 0.5, 0.9)
  expect_equal(qrevweibull(log(p), -1, 2, 0.8, log.p = TRUE),
               qrevweibull(p, -1, 2, 0.8))
  expect_equal(qrevweibull(p, -1, 2, 0.8, lower.tail = FALSE),
               qrevweibull(1 - p, -1, 2, 0.8))
})

test_that("prevweibull: log.p returns the log of the CDF", {
  q <- c(-3, -2, -1)
  expect_equal(prevweibull(q, -1, 0.5, 0.8, log.p = TRUE),
               log(prevweibull(q, -1, 0.5, 0.8)))
})

test_that("qrevweibull and rrevweibull reject a zero scale as d and p do", {
  expect_error(qrevweibull(0.5, scale = 0))
  expect_error(rrevweibull(5, scale = 0))
})

test_that("dnweibull and friends are synonyms of the revweibull functions", {
  x <- c(-3, -2, -1)
  expect_equal(dnweibull(x, -1, 2, 0.8), drevweibull(x, -1, 2, 0.8))
  expect_equal(pnweibull(x, -1, 2, 0.8), prevweibull(x, -1, 2, 0.8))
  expect_equal(qnweibull(0.4, -1, 2, 0.8), qrevweibull(0.4, -1, 2, 0.8))
})

test_that("drevweibull: integrates to 1 (shape=2)", {
  # shape=2, loc=0, scale=1: support (-inf, 0)
  x <- seq(-20, 0, length.out = 100001)
  dx <- x[2] - x[1]
  expect_equal(sum(drevweibull(x, 0, 1, 2)) * dx, 1, tolerance = 1e-3)
})
