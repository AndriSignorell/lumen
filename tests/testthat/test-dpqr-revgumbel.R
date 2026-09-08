library(testthat)
library(lumen)

tol <- 1e-6

test_that("drevgumbel: density >= 0", {
  expect_true(all(drevgumbel(seq(-3, 3, by = 0.25)) >= 0))
})

test_that("drevgumbel: integrates to 1", {
  x <- seq(-20, 10, length.out = 100001)
  dx <- x[2] - x[1]
  expect_equal(sum(drevgumbel(x, loc = 0, scale = 1)) * dx, 1,
               tolerance = 1e-4)
})

test_that("drevgumbel: mode at location", {
  # mode of revgumbel = location; d/dx = 0 at x = location + scale * log(1) = location
  # density at mode = exp(0)*exp(-1)/scale = 1/(e*scale)
  expect_equal(drevgumbel(0, loc = 0, scale = 1), exp(-1), tolerance = tol)
})

test_that("drevgumbel: invalid scale throws error", {
  expect_error(drevgumbel(1, scale = -1))
  expect_error(drevgumbel(1, scale = 0))
})

test_that("prevgumbel: in [0,1]", {
  q <- seq(-5, 5, by = 0.5)
  p <- prevgumbel(q)
  expect_true(all(p >= 0 & p <= 1))
})

test_that("prevgumbel: non-decreasing", {
  q <- seq(-5, 5, by = 0.25)
  expect_true(all(diff(prevgumbel(q)) >= 0))
})

test_that("prevgumbel: relation to pgumbel", {
  # revgumbel(loc, scale) = -Gumbel(-loc, scale)
  # prevgumbel(q) = 1 - pgumbel(-q, loc=-loc, scale=scale)
  q <- c(-2, -1, 0, 1)
  expect_equal(prevgumbel(q, loc = 0, scale = 1),
               1 - pgumbel(-q, loc = 0, scale = 1), tolerance = tol)
})

test_that("qrevgumbel: prevgumbel(qrevgumbel(p)) == p roundtrip", {
  # qrevgumbel(p) = loc + scale*log(-log(1-p)) inverts the CDF
  # prevgumbel(q) = 1 - exp(-exp((q-loc)/scale))
  p <- c(0.1, 0.25, 0.5, 0.75, 0.9)
  expect_equal(prevgumbel(qrevgumbel(p, loc = 1, scale = 2),
                          loc = 1, scale = 2), p, tolerance = tol)
})

test_that("qrevgumbel: median anchor and monotonicity", {
  # pins the orientation of the quantile function,
  # guards against re-inverting it (the median is invariant)
  expect_equal(qrevgumbel(0.5), log(log(2)), tolerance = tol)
  expect_true(all(diff(qrevgumbel(seq(0.05, 0.95, by = 0.05))) > 0))
})

test_that("rrevgumbel: returns correct length", {
  set.seed(1)
  expect_length(rrevgumbel(50), 50)
})


test_that("prevgumbel invalid scale throws error", {
  expect_error(prevgumbel(1, scale = 0))
})

test_that("qrevgumbel invalid scale throws error", {
  expect_error(qrevgumbel(0.5, scale = 0))
})

test_that("qrevgumbelExp equals exp(qrevgumbel())", {
  p <- c(0.2, 0.5, 0.8)
  
  expect_equal(
    qrevgumbelExp(p),
    exp(qrevgumbel(p))
  )
})

test_that("rrevgumbel invalid scale throws error", {
  expect_error(rrevgumbel(10, scale = -1))
})

test_that("drevgumbel: log=TRUE returns the log density", {
  x <- c(-1, 0, 1)
  expect_equal(drevgumbel(x, 1, 2, log = TRUE), log(drevgumbel(x, 1, 2)))
})

test_that("prevgumbel: lower.tail and log.p", {
  q <- c(-1, 0, 1)
  expect_equal(prevgumbel(q, lower.tail = FALSE), 1 - prevgumbel(q))
  expect_equal(prevgumbel(q, log.p = TRUE), log(prevgumbel(q)))
  expect_equal(prevgumbel(q, lower.tail = FALSE, log.p = TRUE),
               log(1 - prevgumbel(q)))
})

test_that("qrevgumbel: p = 0 and p = 1 give the end points of the support", {
  expect_equal(qrevgumbel(c(0, 1)), c(-Inf, Inf))
})

test_that("qrevgumbel: log.p and lower.tail", {
  p <- c(0.1, 0.5, 0.9)
  expect_equal(qrevgumbel(log(p), 1, 2, log.p = TRUE), qrevgumbel(p, 1, 2))
  expect_equal(qrevgumbel(p, 1, 2, lower.tail = FALSE), qrevgumbel(1 - p, 1, 2))
})

test_that("qrevgumbelExp takes the same parameters as qrevgumbel", {
  p <- c(0.2, 0.5, 0.8)
  expect_equal(qrevgumbelExp(p, 1, 2), exp(qrevgumbel(p, 1, 2)))
})

test_that("revgumbel is the reflected Gumbel", {
  # X = loc - scale * Y for a standard Gumbel Y
  x <- c(-2, -1, 0, 1)
  expect_equal(drevgumbel(x, 1, 2), dgumbel(2 - x, 1, 2))
})

