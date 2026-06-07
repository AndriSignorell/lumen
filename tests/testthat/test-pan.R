library(testthat)
library(lumen)

# pan(A, M, C, N) computes P(sum_i lambda_i * chi^2_1 <= x)
# A[1]    = x  (the quantile value, 0-based R index)
# A[2..M+1] = eigenvalues lambda_1..lambda_M  (1-based: A[2..M+1])
# BUT in C++ 0-based: A[0]=x, A[1..M]=eigenvalues, so R passes length M+1 vector
# M  = number of eigenvalues
# C  = coefficient (typically 0 for pure quadratic form)
# N  = number of quadrature points

# Helper: build A vector as expected by pan()
# A[1] = x, A[2:(M+1)] = eigenvalues  (R 1-based)
make_A <- function(x, eigenvalues) c(x, eigenvalues)


test_that("pan: returns value in [0, 1]", {
  eigs <- c(3, 2, 1)
  res  <- pan(make_A(3.0, eigs), M = 3, C = 0, N = 20)
  expect_gte(res, 0)
  expect_lte(res, 1)
})

test_that("pan: very large x gives p close to 1", {
  eigs <- c(3, 2, 1)
  res  <- pan(make_A(1e6, eigs), M = 3, C = 0, N = 20)
  expect_gt(res, 0.99)
})

test_that("pan: very small (negative) x gives p close to 0", {
  eigs <- c(3, 2, 1)
  res  <- pan(make_A(-1e6, eigs), M = 3, C = 0, N = 20)
  expect_lt(res, 0.01)
})

test_that("pan: monotone in x", {
  eigs <- c(3, 2, 1)
  xs   <- c(0, 1, 3, 6, 10)
  ps   <- sapply(xs, function(x) pan(make_A(x, eigs), M = 3, C = 0, N = 20))
  expect_true(all(diff(ps) >= 0))
})

test_that("pan: increasing eigenvalues give monotone CDF", {
  eigs <- c(1.0, 2.0, 3.0)
  xs   <- c(1, 3, 6, 10, 20)
  ps   <- sapply(xs, function(x) pan(make_A(x, eigs), M = 3, C = 0, N = 30))
  expect_true(all(diff(ps) >= 0))
})

test_that("pan: larger sum of eigenvalues shifts CDF right", {
  # same x, bigger eigenvalues -> smaller CDF value
  eigs_small <- c(1.0, 1.0, 1.0)
  eigs_large <- c(3.0, 3.0, 3.0)
  q <- 1.5   # small enough that large eigenvalues haven't saturated CDF
  r_small <- pan(make_A(q, eigs_small), M = 3, C = 0, N = 30)
  r_large <- pan(make_A(q, eigs_large), M = 3, C = 0, N = 30)
  expect_gt(r_small, r_large)
})

test_that("pan: N=1 still returns value in [0,1]", {
  eigs <- c(2, 1)
  res  <- pan(make_A(2.0, eigs), M = 2, C = 0, N = 1)
  expect_gte(res, 0)
  expect_lte(res, 1)
})

test_that("pan: larger N gives more accurate result", {
  # for equal eigenvalues, compare against pchisq
  k    <- 3; lam <- 1.5; q <- 6
  ref  <- pchisq(q / lam, df = k)
  eigs <- rep(lam, k)
  r10  <- abs(pan(make_A(q, eigs), M = k, C = 0, N = 10)  - ref)
  r100 <- abs(pan(make_A(q, eigs), M = k, C = 0, N = 100) - ref)
  expect_lte(r100, r10 + 1e-10)   # larger N not worse
})

test_that("pan: matches lmtest::dwtest p-value (via durbinWatsonTest)", {
  skip_if_not_installed("lmtest")
  set.seed(1)
  n  <- 25
  x  <- rnorm(n)
  y  <- 1 + x + rnorm(n)
  fit <- lm(y ~ x)
  
  ref <- lmtest::dwtest(fit, alternative = "greater")$p.value
  res <- durbinWatsonTest(fit, alternative = "greater", exact = TRUE)$p.value
  expect_equal(res, ref, tolerance = 1e-4)
})
