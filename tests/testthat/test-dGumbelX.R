library(testthat)
library(lumen)

tol <- 1e-6

# Parameters: loc1 <= loc2 required
L1 <- 0; S1 <- 1.1; L2 <- 1; S2 <- 0.5

test_that("dgumbelx: density >= 0", {
  expect_true(all(dgumbelx(seq(-1, 6, by = 0.5), L1, S1, L2, S2) >= 0))
})

test_that("dgumbelx: integrates to 1", {
  x <- seq(-5, 15, length.out = 100001)
  dx <- x[2] - x[1]
  expect_equal(sum(dgumbelx(x, L1, S1, L2, S2)) * dx, 1, tolerance = 1e-3)
})

test_that("dgumbelx: log=TRUE", {
  x <- c(1, 2, 3)
  expect_equal(dgumbelx(x, L1, S1, L2, S2, log = TRUE),
               log(dgumbelx(x, L1, S1, L2, S2)), tolerance = tol)
})

test_that("dgumbelx: loc1 > loc2 throws error", {
  expect_error(dgumbelx(1, loc1 = 2, loc2 = 1))
})

test_that("pgumbelx: in [0,1] and non-decreasing", {
  q <- seq(-2, 10, by = 0.5)
  p <- pgumbelx(q, L1, S1, L2, S2)
  expect_true(all(p >= 0 & p <= 1))
  expect_true(all(diff(p) >= 0))
})

test_that("pgumbelx: lower.tail=FALSE complement", {
  q <- c(0, 1, 3)
  expect_equal(pgumbelx(q, L1, S1, L2, S2, lower.tail = FALSE),
               1 - pgumbelx(q, L1, S1, L2, S2), tolerance = tol)
})

test_that("pgumbelx: when loc1=loc2 and scale1=scale2, equals pgumbel^2", {
  # max of two identical Gumbels: CDF = F(q)^2
  q <- c(0, 1, 2)
  expect_equal(pgumbelx(q, 0, 1, 0, 1),
               pgumbel(q)^2, tolerance = tol)
})

test_that("pgumbelx(qgumbelx(p)) roundtrip", {
  p <- c(0.2, 0.5, 0.8)
  q <- qgumbelx(p, interval = c(-5, 20), L1, S1, L2, S2)
  expect_equal(pgumbelx(q, L1, S1, L2, S2), p, tolerance = 1e-5)
})

test_that("rgumbelx: returns correct length", {
  set.seed(1)
  expect_length(rgumbelx(50, L1, S1, L2, S2), 50)
})

test_that("rgumbelx: values >= max(loc1, loc2) not guaranteed but finite", {
  set.seed(1)
  r <- rgumbelx(100, L1, S1, L2, S2)
  expect_true(all(is.finite(r)))
})
