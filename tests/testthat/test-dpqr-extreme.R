library(testthat)
library(lumen)

tol <- 1e-6

# --- dextreme (max of mlen iid rvs) ---

test_that("dextreme mlen=1 equals base density", {
  x <- c(1, 2, 3)
  expect_equal(dextreme(x, distn = "norm", mlen = 1),
               dnorm(x), tolerance = tol)
})

test_that("dextreme: density >= 0", {
  expect_true(all(dextreme(seq(-3, 3, by = 0.5), distn = "norm", mlen = 3) >= 0))
})

test_that("dextreme: integrates to 1", {
  x <- seq(-5, 10, length.out = 100001)
  dx <- x[2] - x[1]
  expect_equal(sum(dextreme(x, distn = "norm", mean = 0.5, sd = 1.2, mlen = 5)) * dx,
               1, tolerance = 1e-3)
})

test_that("dextreme log=TRUE", {
  x <- c(1, 2, 3)
  expect_equal(dextreme(x, distn = "norm", mlen = 3, log = TRUE),
               log(dextreme(x, distn = "norm", mlen = 3)), tolerance = tol)
})

test_that("dextreme largest=FALSE (minima): integrates to 1", {
  x <- seq(-5, 10, length.out = 100001)
  dx <- x[2] - x[1]
  expect_equal(sum(dextreme(x, distn = "exp", mlen = 2, largest = FALSE)) * dx,
               1, tolerance = 1e-3)
})

test_that("dextreme: invalid mlen throws error", {
  expect_error(dextreme(1, distn = "norm", mlen = 0))
  expect_error(dextreme(1, distn = "norm", mlen = 1.5))
})

# --- pextreme ---

test_that("pextreme mlen=1 equals base CDF", {
  q <- c(-1, 0, 1)
  expect_equal(pextreme(q, distn = "norm", mlen = 1), pnorm(q), tolerance = tol)
})

test_that("pextreme: CDF = F(q)^mlen", {
  q <- c(0, 1, 2)
  expect_equal(pextreme(q, distn = "norm", mlen = 3),
               pnorm(q)^3, tolerance = tol)
})

test_that("pextreme: non-decreasing", {
  q <- seq(-3, 3, by = 0.25)
  expect_true(all(diff(pextreme(q, distn = "norm", mlen = 5)) >= 0))
})

test_that("pextreme: lower.tail=FALSE complement", {
  q <- c(0, 1, 2)
  expect_equal(pextreme(q, distn = "exp", rate = 1.2, mlen = 2, lower.tail = FALSE),
               1 - pextreme(q, distn = "exp", rate = 1.2, mlen = 2), tolerance = tol)
})

# --- qextreme ---

test_that("pextreme(qextreme(p)) roundtrip", {
  p <- c(0.1, 0.5, 0.9)
  expect_equal(
    pextreme(qextreme(p, distn = "exp", rate = 1.2, mlen = 2),
             distn = "exp", rate = 1.2, mlen = 2),
    p, tolerance = tol)
})

test_that("qextreme: p = 0 and p = 1 give the end points of the support", {
  expect_equal(qextreme(c(0, 1), distn = "norm", mlen = 2), c(-Inf, Inf))
})

test_that("qextreme: p outside [0,1] gives NaN with a warning", {
  expect_warning(res <- qextreme(1.1, distn = "norm", mlen = 2), "NaN")
  expect_true(is.nan(res))
})

test_that("qextreme: log.p and lower.tail", {
  p <- c(0.1, 0.5, 0.9)
  expect_equal(qextreme(log(p), distn = "exp", rate = 1.2, mlen = 2,
                        log.p = TRUE),
               qextreme(p, distn = "exp", rate = 1.2, mlen = 2))
  expect_equal(qextreme(p, distn = "exp", rate = 1.2, mlen = 2,
                        lower.tail = FALSE),
               qextreme(1 - p, distn = "exp", rate = 1.2, mlen = 2))
})

test_that("pextreme: log.p returns the log of the CDF", {
  q <- c(0, 1, 2)
  expect_equal(pextreme(q, distn = "norm", mlen = 3, log.p = TRUE),
               log(pextreme(q, distn = "norm", mlen = 3)))
})

test_that("dextreme mlen = 1: finite where the CDF underflows to zero", {
  # (mlen - 1) * log(F) once produced 0 * -Inf = NaN here
  expect_equal(dextreme(-40, distn = "norm", mlen = 1), dnorm(-40))
  expect_equal(dextreme(0, distn = "unif", mlen = 1), dunif(0))
  expect_false(is.nan(dextreme(-40, distn = "norm", mlen = 1, log = TRUE)))
})
