library(testthat)
library(lumen)

tol <- 1e-6

# --- ddirichlet ---

test_that("ddirichlet: symmetric Dirichlet(1,1,1) = 2 on simplex", {
  # Dir(1,1,1) is uniform on simplex; density = Gamma(3)/Gamma(1)^3 = 2
  x <- matrix(c(0.2, 0.3, 0.5), nrow = 1)
  expect_equal(ddirichlet(x, alpha = c(1,1,1)), 2, tolerance = tol)
})

test_that("ddirichlet: density >= 0", {
  set.seed(1)
  x <- rdirichlet(20, c(2, 3, 1))
  d <- ddirichlet(x, c(2, 3, 1))
  expect_true(all(d >= 0))
})

test_that("ddirichlet: vector input (single draw)", {
  x <- c(0.2, 0.3, 0.5)
  d <- ddirichlet(x, c(1, 1, 1))
  expect_equal(d, 2, tolerance = tol)
})

test_that("ddirichlet: off-simplex gives 0", {
  x <- matrix(c(0.2, 0.3, 0.6), nrow = 1)  # sums to 1.1
  expect_equal(ddirichlet(x, c(1, 1, 1)), 0)
})

test_that("ddirichlet: negative values give 0", {
  x <- matrix(c(-0.1, 0.6, 0.5), nrow = 1)
  expect_equal(ddirichlet(x, c(1, 1, 1)), 0)
})

test_that("ddirichlet log=TRUE", {
  x <- c(0.2, 0.3, 0.5)
  expect_equal(ddirichlet(x, c(2, 3, 4), log = TRUE),
               log(ddirichlet(x, c(2, 3, 4))), tolerance = tol)
})

test_that("ddirichlet: alpha <= 0 throws error", {
  expect_error(ddirichlet(c(0.2, 0.3, 0.5), c(1, 0, 1)))
  expect_error(ddirichlet(c(0.2, 0.3, 0.5), c(1, -1, 1)))
})

test_that("ddirichlet: mismatched length throws error", {
  expect_error(ddirichlet(c(0.2, 0.3, 0.5), c(1, 1)))
})

# --- rdirichlet ---

test_that("rdirichlet: rows sum to 1", {
  set.seed(42)
  x <- rdirichlet(100, c(1, 2, 3))
  expect_equal(rowSums(x), rep(1, 100), tolerance = 1e-10)
})

test_that("rdirichlet: all values in [0,1]", {
  set.seed(1)
  x <- rdirichlet(50, c(1, 1, 1, 1))
  expect_true(all(x >= 0 & x <= 1))
})

test_that("rdirichlet: returns n rows and k columns", {
  x <- rdirichlet(30, c(1, 2, 3))
  expect_equal(nrow(x), 30)
  expect_equal(ncol(x), 3)
})

test_that("rdirichlet: alpha <= 0 throws error", {
  expect_error(rdirichlet(10, c(1, 0, 1)))
})

# --- pdirichlet ---

test_that("pdirichlet: probability in [0,1]", {
  set.seed(1)
  p <- pdirichlet(c(0.2, 0.3, 0.5), c(1, 1, 1), n_sim = 1e4)
  expect_true(p >= 0 && p <= 1)
})

test_that("pdirichlet: P(X <= 1) = 1 (maximal point on simplex)", {
  p <- pdirichlet(c(1, 1, 1), c(1, 1, 1), n_sim = 1e5)
  expect_equal(p, 1, tolerance = 0.01)
})

# --- qdirichlet ---

test_that("qdirichlet: always throws error (not defined)", {
  expect_error(qdirichlet())
})
