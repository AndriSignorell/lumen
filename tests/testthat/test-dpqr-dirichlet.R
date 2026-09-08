library(testthat)
library(lumen)

tol <- 1e-6

# --- ddirichlet ---

test_that("ddirichlet: symmetric Dirichlet(1,1,1) = 2 on simplex", {
  # Dir(1,1,1) is uniform on simplex; density = Gamma(3)/Gamma(1)^3 = 2
  x <- matrix(c(0.2, 0.3, 0.5), nrow = 1)
  expect_equal(ddirichlet(x, concentration = c(1,1,1)), 2, tolerance = tol)
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

test_that("ddirichlet: concentration <= 0 throws error", {
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

test_that("rdirichlet: concentration <= 0 throws error", {
  expect_error(rdirichlet(10, c(1, 0, 1)))
})

# --- pdirichlet ---

test_that("pdirichlet: probability in [0,1]", {
  set.seed(1)
  p <- pdirichlet(c(0.5, 0.6, 0.7), c(1, 1, 1), R = 1e4)
  expect_true(p >= 0 && p <= 1)
})

test_that("pdirichlet: agrees with an exactly known value", {
  # for Dir(1,1,1), P(X_i > q_i) = (1 - q_i)^2, and the three events are
  # disjoint once the thresholds sum above 1, so P = 1 - sum (1 - q_i)^2
  q <- c(0.5, 0.6, 0.7)
  set.seed(1)
  expect_equal(pdirichlet(q, c(1, 1, 1), R = 2e5), 1 - sum((1 - q)^2),
               tolerance = 0.01)
})

test_that("pdirichlet: a degenerate region has probability zero", {
  # x1 <= 0.2 and x2 <= 0.3 force x1 + x2 <= 0.5, while x3 <= 0.5 forces
  # x1 + x2 >= 0.5: the region is a face of the simplex
  set.seed(1)
  expect_equal(pdirichlet(c(0.2, 0.3, 0.5), c(1, 1, 1), R = 1e4), 0)
})

test_that("pdirichlet: set.seed makes the simulation reproducible", {
  q <- c(0.5, 0.6, 0.7)
  set.seed(42); a <- pdirichlet(q, c(1, 1, 1), R = 5e4)
  set.seed(42); b <- pdirichlet(q, c(1, 1, 1), R = 5e4)
  set.seed(43); d <- pdirichlet(q, c(1, 1, 1), R = 5e4)
  expect_identical(a, b)
  expect_false(identical(a, d))
})

test_that("pdirichlet: P(X <= 1) = 1 (maximal point on simplex)", {
  p <- pdirichlet(c(1, 1, 1), c(1, 1, 1), R = 1e5)
  expect_equal(p, 1, tolerance = 0.01)
})

# --- qdirichlet ---

test_that("qdirichlet: always throws error (not defined)", {
  expect_error(qdirichlet())
  # the arguments are accepted so that the message is reached
  expect_error(qdirichlet(0.5, c(1, 1, 1)), "no quantile function")
})

test_that("pdirichlet: mismatched length throws error", {
  expect_error(pdirichlet(c(0.2, 0.8), c(1, 1, 1), R = 1000))
})

test_that("pdirichlet: invalid concentration or R throws error", {
  expect_error(pdirichlet(c(0.2, 0.3, 0.5), c(1, 0, 1)))
  expect_error(pdirichlet(c(0.2, 0.3, 0.5), c(1, 1, 1), R = 0))
})
