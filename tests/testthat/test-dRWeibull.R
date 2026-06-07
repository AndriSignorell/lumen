library(testthat)
library(lumen)

tol <- 1e-6

test_that("drweibull: density >= 0", {
  expect_true(all(drweibull(seq(-5, -0.01, by = 0.1), 0, 1, 1) >= 0))
})

test_that("drweibull: density = 0 for x >= loc", {
  expect_equal(drweibull(c(0, 1, 2), loc = 0), c(0, 0, 0))
})

test_that("drweibull: integrates to 1", {
  x <- seq(-30, 0, length.out = 100001)
  dx <- x[2] - x[1]
  expect_equal(sum(drweibull(x, 0, 1, 1)) * dx, 1, tolerance = 1e-3)
})

test_that("drweibull log=TRUE", {
  x <- c(-3, -2, -1)
  expect_equal(drweibull(x, -1, 0.5, 0.8, log = TRUE),
               log(drweibull(x, -1, 0.5, 0.8)), tolerance = tol)
})

test_that("drweibull: invalid scale/shape throws error", {
  expect_error(drweibull(-1, scale = 0))
  expect_error(drweibull(-1, shape = 0))
})

test_that("prweibull: in [0,1] and non-decreasing", {
  q <- seq(-10, 0, by = 0.5)
  p <- prweibull(q, 0, 1, 1)
  expect_true(all(p >= 0 & p <= 1))
  expect_true(all(diff(p) >= 0))
})

test_that("prweibull: lower.tail=FALSE complement", {
  q <- c(-3, -2, -1)
  expect_equal(prweibull(q, -1, 0.5, 0.8, lower.tail = FALSE),
               1 - prweibull(q, -1, 0.5, 0.8), tolerance = tol)
})

test_that("prweibull(qrweibull(p)) roundtrip", {
  p <- c(0.1, 0.25, 0.5, 0.75, 0.9)
  expect_equal(prweibull(qrweibull(p, -1, 2, 0.8), -1, 2, 0.8), p, tolerance = tol)
})

test_that("qrweibull: invalid p throws error", {
  expect_error(qrweibull(0))
  expect_error(qrweibull(1))
})

test_that("drweibull: integrates to 1 (shape=2)", {
  # shape=2, loc=0, scale=1: support (-inf, 0)
  x <- seq(-20, 0, length.out = 100001)
  dx <- x[2] - x[1]
  expect_equal(sum(drweibull(x, 0, 1, 2)) * dx, 1, tolerance = 1e-3)
})
