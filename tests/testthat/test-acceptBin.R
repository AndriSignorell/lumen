library(testthat)
library(lumen)

# acceptBin(x, n, p) returns the "acceptability" of p for the Blaker
# exact CI. It is always in (0, 1], and p is *outside* the CI when
# acceptBin(x, n, p) < alpha.

test_that("acceptBin: returns value in (0, 1]", {
  res <- acceptBin(5L, 20L, 0.25)
  expect_gt(res, 0)
  expect_lte(res, 1)
})

test_that("acceptBin: p = x/n (MLE) gives large acceptability", {
  # at the MLE the observed value is central; acceptability should be large
  res <- acceptBin(10L, 20L, 0.5)
  expect_gt(res, 0.5)
})

test_that("acceptBin: extreme p gives small acceptability", {
  # p = 0.99 implausible for x=2, n=20
  res <- acceptBin(2L, 20L, 0.99)
  expect_lt(res, 0.05)
})

test_that("acceptBin: x = 0 and p close to 0 gives large acceptability", {
  res <- acceptBin(0L, 20L, 0.01)
  expect_gt(res, 0.05)
})

test_that("acceptBin: x = n and p close to 1 gives large acceptability", {
  res <- acceptBin(20L, 20L, 0.99)
  expect_gt(res, 0.05)
})

test_that("acceptBin: x = 0 and p = 0.5 gives small acceptability", {
  res <- acceptBin(0L, 20L, 0.5)
  expect_lt(res, 0.05)
})

test_that("acceptBin: symmetry - acceptBin(x,n,p) == acceptBin(n-x,n,1-p)", {
  res1 <- acceptBin(3L,  15L, 0.3)
  res2 <- acceptBin(12L, 15L, 0.7)
  expect_equal(res1, res2, tolerance = 1e-12)
})

test_that("acceptBin: result is numeric scalar", {
  res <- acceptBin(5L, 10L, 0.5)
  expect_true(is.numeric(res))
  expect_length(res, 1L)
})
