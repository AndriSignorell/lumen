library(testthat)
library(lumen)

tol <- 1e-6

test_that("dfrechet: density >= 0", {
  expect_true(all(dfrechet(seq(0.1, 5, by = 0.1), 0, 1, 1) >= 0))
})

test_that("dfrechet: density = 0 for x <= loc", {
  expect_equal(dfrechet(c(-1, 0, 1), loc = 1, scale = 1, shape = 1),
               c(0, 0, 0))
})

test_that("dfrechet: integrates to 1", {
  # shape=2: lighter tail, upper bound 200 is sufficient
  x <- seq(0.001, 200, length.out = 500001)
  dx <- x[2] - x[1]
  expect_equal(sum(dfrechet(x, 0, 1, 2)) * dx, 1, tolerance = 1e-3)
})

test_that("dfrechet log=TRUE", {
  x <- c(1.5, 2, 3)
  expect_equal(dfrechet(x, 1, 0.5, 0.8, log = TRUE),
               log(dfrechet(x, 1, 0.5, 0.8)), tolerance = tol)
})

test_that("dfrechet: invalid scale/shape throws error", {
  expect_error(dfrechet(2, scale = 0))
  expect_error(dfrechet(2, shape = 0))
})

test_that("pfrechet: in [0,1] and non-decreasing", {
  q <- seq(0.1, 10, by = 0.5)
  p <- pfrechet(q, 0, 1, 1)
  expect_true(all(p >= 0 & p <= 1))
  expect_true(all(diff(p) >= 0))
})

test_that("pfrechet: lower.tail=FALSE complement", {
  q <- c(1, 2, 3)
  expect_equal(pfrechet(q, 1, 0.5, 0.8, lower.tail = FALSE),
               1 - pfrechet(q, 1, 0.5, 0.8), tolerance = tol)
})

test_that("pfrechet(qfrechet(p)) roundtrip", {
  p <- c(0.1, 0.25, 0.5, 0.75, 0.9)
  expect_equal(pfrechet(qfrechet(p, 1, 2, 0.8), 1, 2, 0.8), p, tolerance = tol)
})

test_that("qfrechet: invalid p throws error", {
  expect_error(qfrechet(0))
  expect_error(qfrechet(1))
})
