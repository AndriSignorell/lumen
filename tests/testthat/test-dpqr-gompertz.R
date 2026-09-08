library(testthat)
library(lumen)

tol <- 1e-6

test_that("dgompertz: density >= 0", {
  expect_true(all(dgompertz(seq(0, 5, by = 0.25), shape = 1, rate = 1) >= 0))
})

test_that("dgompertz shape=0: equals dexp(x, rate=rate)", {
  x <- c(0.5, 1, 2, 3)
  expect_equal(dgompertz(x, shape = 0, rate = 1),
               dexp(x, rate = 1), tolerance = tol)
})

test_that("dgompertz: integrates to 1 (shape > 0)", {
  x <- seq(0, 20, length.out = 100001)
  dx <- x[2] - x[1]
  expect_equal(sum(dgompertz(x, shape = 1, rate = 1)) * dx, 1, tolerance = 1e-3)
})

test_that("dgompertz log=TRUE", {
  x <- c(0.5, 1, 2)
  expect_equal(dgompertz(x, 1, 1, log = TRUE),
               log(dgompertz(x, 1, 1)), tolerance = tol)
})

test_that("dgompertz: negative x gives 0", {
  expect_equal(dgompertz(-1, shape = 1, rate = 1), 0)
})

test_that("pgompertz shape=0: equals pexp", {
  q <- c(0.5, 1, 2, 3)
  expect_equal(pgompertz(q, shape = 0, rate = 1),
               pexp(q, rate = 1), tolerance = tol)
})

test_that("pgompertz: in [0,1] and non-decreasing", {
  q <- seq(0, 5, by = 0.25)
  p <- pgompertz(q, shape = 1, rate = 1)
  expect_true(all(p >= 0 & p <= 1))
  expect_true(all(diff(p) >= 0))
})

test_that("pgompertz: lower.tail=FALSE complement", {
  q <- c(0.5, 1, 2)
  expect_equal(pgompertz(q, 1, 1, lower.tail = FALSE),
               1 - pgompertz(q, 1, 1), tolerance = tol)
})

test_that("dgompertz: an infinite argument gives a zero density", {
  expect_equal(dgompertz(Inf, shape = 1, rate = 1), 0)
})

test_that("dgompertz: overflow of exp(shape * x) does not give NaN", {
  expect_equal(dgompertz(1e6, shape = 1, rate = 1), 0)
  expect_equal(pgompertz(1e6, shape = 1, rate = 1), 1)
})

test_that("gompertz: NA and NaN are kept apart", {
  expect_true(is.na(dgompertz(NA_real_, 1, 1)) &&
              !is.nan(dgompertz(NA_real_, 1, 1)))
  expect_true(is.nan(dgompertz(NaN, 1, 1)))
  expect_true(is.nan(pgompertz(NaN, 1, 1)))
})

test_that("gompertz: a non-positive rate gives NaN with one warning", {
  expect_warning(d <- dgompertz(1:5, 1, rep(-1, 5)), "rate")
  expect_true(all(is.nan(d)))
  expect_warning(p <- pgompertz(1, 1, 0), "rate")
  expect_true(is.nan(p))
  # d, p, q and r agree on what an admissible parameter is
  expect_warning(q <- qgompertz(0.5, 1, 0), "rate")
  expect_true(is.nan(q))
})

test_that("pgompertz: a negative shape leaves mass at infinity", {
  # P(X = Inf) = exp(rate/shape), so the CDF tends to 1 - exp(rate/shape)
  expect_equal(pgompertz(Inf, shape = -0.5, rate = 1), 1 - exp(-2),
               tolerance = 1e-10)
  expect_equal(pgompertz(1e3, shape = -0.5, rate = 1),
               pgompertz(Inf, shape = -0.5, rate = 1), tolerance = 1e-10)
  expect_equal(qgompertz(0.9, shape = -0.5, rate = 1), Inf)
})

test_that("gompertz: a shape below the tolerance is the exponential limit", {
  expect_equal(pgompertz(2, shape = 1e-13, rate = 1), pexp(2, 1))
  expect_equal(qgompertz(0.4, shape = 1e-13, rate = 1), qexp(0.4, 1))
})

test_that("pgompertz(qgompertz(p)) roundtrip", {
  p <- c(0.1, 0.25, 0.5, 0.75, 0.9)
  expect_equal(pgompertz(qgompertz(p, shape = 1, rate = 1), shape = 1, rate = 1),
               p, tolerance = tol)
})
