library(testthat)
library(lumen)

tol <- 1e-6

test_that("dorder j=mlen=1: equals base density", {
  x <- c(0, 1, 2)
  expect_equal(dorder(x, distn = "norm", mlen = 1, j = 1),
               dnorm(x), tolerance = tol)
})

test_that("dorder: density >= 0", {
  expect_true(all(dorder(seq(-3, 3, by = 0.5), distn = "norm",
                         mlen = 5, j = 2) >= 0))
})

test_that("dorder: integrates to 1", {
  x <- seq(-6, 6, length.out = 100001)
  dx <- x[2] - x[1]
  expect_equal(sum(dorder(x, distn = "norm", mlen = 5, j = 2)) * dx,
               1, tolerance = 1e-3)
})

test_that("dorder log=TRUE", {
  x <- c(1, 2, 3)
  expect_equal(dorder(x, distn = "norm", mlen = 5, j = 2, log = TRUE),
               log(dorder(x, distn = "norm", mlen = 5, j = 2)), tolerance = tol)
})

test_that("dorder: j > mlen throws error", {
  expect_error(dorder(1, distn = "norm", mlen = 3, j = 4))
})

test_that("dorder: invalid mlen throws error", {
  expect_error(dorder(1, distn = "norm", mlen = 0))
})

# --- porder ---

test_that("porder j=mlen=1: equals base CDF", {
  q <- c(-1, 0, 1)
  expect_equal(porder(q, distn = "norm", mlen = 1, j = 1),
               pnorm(q), tolerance = tol)
})

test_that("porder: in [0,1] and non-decreasing", {
  q <- seq(-3, 3, by = 0.25)
  p <- porder(q, distn = "norm", mlen = 5, j = 3)
  expect_true(all(p >= 0 & p <= 1))
  expect_true(all(diff(p) >= 0))
})

test_that("porder largest=FALSE: CDF of minimum", {
  # P(X_{(1)} <= q) = 1 - (1-F(q))^mlen
  q <- c(-1, 0, 1)
  expect_equal(porder(q, distn = "norm", mlen = 3, j = 1, largest = FALSE),
               1 - (1 - pnorm(q))^3, tolerance = tol)
})

# --- rorder ---

test_that("rorder: returns correct length", {
  set.seed(1)
  expect_length(rorder(20, distn = "norm", mlen = 5, j = 2), 20)
})

test_that("rorder: values within plausible range", {
  set.seed(1)
  r <- rorder(500, distn = "norm", mlen = 5, j = 3)
  expect_true(all(is.finite(r)))
  # median of 3rd-largest of 5 standard normals should be around 0
  expect_equal(median(r), 0, tolerance = 0.1)
})
