library(testthat)
library(lumen)

tol <- 1e-6

test_that("dgpd: density >= 0", {
  expect_true(all(dgpd(seq(0.01, 5, by = 0.1), 0, 1, 0.5) >= 0))
})

test_that("dgpd shape=0: exponential density at loc=0, scale=1", {
  # GPD(shape=0) = Exp(rate=1/scale), density = (1/scale)*exp(-x/scale)
  x <- c(0.5, 1, 2)
  expect_equal(dgpd(x, loc = 0, scale = 1, shape = 0),
               dexp(x, rate = 1), tolerance = tol)
})

test_that("dgpd: integrates to 1 (shape > 0)", {
  # shape=0.3: S(100) = (1+0.3*100)^(-1/0.3) < 1e-4; seq to 100 is sufficient
  x <- seq(0, 100, length.out = 100001)
  dx <- x[2] - x[1]
  expect_equal(sum(dgpd(x, 0, 1, 0.3)) * dx, 1, tolerance = 1e-3)
})

test_that("dgpd log=TRUE", {
  x <- c(2, 3, 4)
  expect_equal(dgpd(x, 1, 0.5, 0.8, log = TRUE),
               log(dgpd(x, 1, 0.5, 0.8)), tolerance = tol)
})

test_that("dgpd: invalid scale throws error", {
  expect_error(dgpd(1, scale = -1))
})

test_that("pgpd shape=0: equals pexp", {
  q <- c(0.5, 1, 2, 3)
  expect_equal(pgpd(q, loc = 0, scale = 1, shape = 0),
               pexp(q, rate = 1), tolerance = tol)
})

test_that("pgpd: in [0,1] and non-decreasing", {
  q <- seq(0, 10, by = 0.5)
  p <- pgpd(q, 0, 1, 0.5)
  expect_true(all(p >= 0 & p <= 1))
  expect_true(all(diff(p) >= 0))
})

test_that("pgpd: lower.tail=FALSE complement", {
  q <- c(1, 2, 3)
  expect_equal(pgpd(q, 1, 0.5, 0.8, lower.tail = FALSE),
               1 - pgpd(q, 1, 0.5, 0.8), tolerance = tol)
})

test_that("pgpd(qgpd(p)) roundtrip", {
  p <- c(0.1, 0.25, 0.5, 0.75, 0.9)
  expect_equal(pgpd(qgpd(p, 1, 2, 0.8), 1, 2, 0.8), p, tolerance = tol)
})

test_that("qgpd: invalid p throws error", {
  expect_error(qgpd(0))
  expect_error(qgpd(1))
})
