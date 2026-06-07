library(testthat)
library(lumen)

test_that("fisherZ: r=0 gives z=0", {
  expect_equal(fisherZ(0), 0)
})

test_that("fisherZ: known value r=0.5", {
  expect_equal(fisherZ(0.5), atanh(0.5), tolerance = 1e-10)
})

test_that("fisherZ: r=1 gives Inf", {
  expect_equal(fisherZ(1), Inf)
})

test_that("fisherZ: r=-1 gives -Inf", {
  expect_equal(fisherZ(-1), -Inf)
})

test_that("fisherZ: antisymmetric", {
  r <- c(0.1, 0.3, 0.7, 0.9)
  expect_equal(fisherZ(-r), -fisherZ(r))
})

test_that("fisherZ vectorises", {
  r <- seq(-0.9, 0.9, by = 0.1)
  expect_length(fisherZ(r), length(r))
})

test_that("fisherZInv: z=0 gives r=0", {
  expect_equal(fisherZInv(0), 0)
})

test_that("fisherZInv: known value", {
  expect_equal(fisherZInv(atanh(0.5)), 0.5, tolerance = 1e-10)
})

test_that("fisherZInv: inverse of fisherZ", {
  r <- seq(-0.9, 0.9, by = 0.1)
  expect_equal(fisherZInv(fisherZ(r)), r, tolerance = 1e-10)
})

test_that("fisherZ: inverse of fisherZInv", {
  z <- seq(-2, 2, by = 0.25)
  expect_equal(fisherZ(fisherZInv(z)), z, tolerance = 1e-10)
})

test_that("fisherZInv: output in [-1, 1]", {
  z <- c(-5, -2, 0, 2, 5)
  r <- fisherZInv(z)
  expect_true(all(r >= -1 & r <= 1))
})
