library(testthat)
library(lumen)

tol <- 1e-6

# --- dgev ---

test_that("dgev shape=0 equals Gumbel density", {
  x <- c(-1, 0, 1, 2)
  expect_equal(dgev(x, loc = 1, scale = 0.5, shape = 0),
               dgumbel(x, loc = 1, scale = 0.5), tolerance = tol)
})

test_that("dgev: density >= 0", {
  expect_true(all(dgev(seq(-5, 10, by = 0.5), shape = 0.5) >= 0))
})

test_that("dgev: integrates to 1 (shape > 0)", {
  # shape=0.5, loc=1, scale=2: lower bound = loc-scale/shape = 1-4 = -3; heavy upper tail
  x <- seq(-3, 500, length.out = 500001)
  dx <- x[2] - x[1]
  expect_equal(sum(dgev(x, loc = 1, scale = 2, shape = 0.5)) * dx, 1, tolerance = 1e-3)
})

test_that("dgev: integrates to 1 (shape < 0)", {
  # shape=-0.5, loc=0, scale=1: support is (-inf, loc-scale/shape) = (-inf, 2)
  x <- seq(-50, 2, length.out = 500001)
  dx <- x[2] - x[1]
  expect_equal(sum(dgev(x, loc = 0, scale = 1, shape = -0.5)) * dx, 1, tolerance = 1e-3)
})

test_that("dgev log=TRUE", {
  x <- c(2, 3, 4)
  expect_equal(dgev(x, 1, 0.5, 0.8, log = TRUE),
               log(dgev(x, 1, 0.5, 0.8)), tolerance = tol)
})

test_that("dgev: invalid scale throws error", {
  expect_error(dgev(1, scale = -1))
})

# --- pgev ---

test_that("pgev: shape=0 equals pgumbel", {
  q <- c(-1, 0, 1, 2)
  expect_equal(pgev(q, loc = 0, scale = 1, shape = 0),
               pgumbel(q), tolerance = tol)
})

test_that("pgev: in [0,1]", {
  expect_true(all(between <- pgev(seq(-3, 5), shape = 0.5) >= 0))
  expect_true(all(pgev(seq(-3, 5), shape = 0.5) <= 1))
})

test_that("pgev: non-decreasing", {
  q <- seq(-3, 5, by = 0.25)
  expect_true(all(diff(pgev(q, shape = 0.8)) >= 0))
})

test_that("pgev: lower.tail=FALSE complement", {
  q <- c(0, 1, 2)
  expect_equal(pgev(q, shape = 0.5, lower.tail = FALSE),
               1 - pgev(q, shape = 0.5), tolerance = tol)
})

# --- qgev ---

test_that("pgev(qgev(p)) roundtrip", {
  p <- c(0.1, 0.25, 0.5, 0.75, 0.9)
  expect_equal(pgev(qgev(p, 1, 2, 0.8), 1, 2, 0.8), p, tolerance = tol)
})

test_that("pgev(qgev(p)) roundtrip shape=0", {
  p <- c(0.1, 0.5, 0.9)
  expect_equal(pgev(qgev(p, 1, 2, 0), 1, 2, 0), p, tolerance = tol)
})

test_that("qgev: invalid p throws error", {
  expect_error(qgev(0))
  expect_error(qgev(1))
})


test_that("dgev shape length > 1 throws error", {
  expect_error(dgev(1, shape = c(0, 1)))
})

test_that("dgev outside support returns zero", {
  expect_equal(
    dgev(-10, loc = 0, scale = 1, shape = 0.5),
    0
  )
})

test_that("pgev shape length > 1 throws error", {
  expect_error(pgev(1, shape = c(0, 1)))
})

test_that("qgev lower.tail=FALSE works", {
  p <- c(0.1, 0.5, 0.9)
  
  expect_equal(
    pgev(qgev(p, lower.tail = FALSE),
         lower.tail = FALSE),
    p,
    tolerance = 1e-6
  )
})

test_that("qgev shape length > 1 throws error", {
  expect_error(qgev(0.5, shape = c(0, 1)))
})

test_that("rgev returns correct length", {
  expect_length(rgev(25), 25)
})

test_that("rgev invalid shape length throws error", {
  expect_error(rgev(10, shape = c(0, 1)))
})

