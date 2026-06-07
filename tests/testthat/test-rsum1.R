library(testthat)
library(lumen)

test_that("rsum1: returns vector of length n", {
  set.seed(1); expect_length(rsum1(5), 5L)
})

test_that("rsum1: values sum to 1", {
  set.seed(1); expect_equal(sum(rsum1(5)), 1, tolerance = 1e-10)
})

test_that("rsum1: all values >= 0", {
  set.seed(1); expect_true(all(rsum1(10) >= 0))
})

test_that("rsum1: n=1 returns 1", {
  set.seed(1); expect_equal(rsum1(1), 1, tolerance = 1e-10)
})

test_that("rsum1: digits rounds and still sums to 1", {
  set.seed(42)
  x <- rsum1(5, digits = 2)
  expect_equal(sum(x), 1, tolerance = 1e-10)
  # all values have at most 2 decimal places
  expect_true(all(abs(x - round(x, 2)) < 1e-10))
})

test_that("rsum1: different seeds give different results", {
  set.seed(1); x1 <- rsum1(5)
  set.seed(2); x2 <- rsum1(5)
  expect_false(isTRUE(all.equal(x1, x2)))
})
