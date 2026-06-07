library(testthat)
library(lumen)

# Analytical values from the docstring examples

test_that("dtri: known values min=10, max=15, mode=12", {
  expect_equal(dtri(12:14, 10, 15, 12),
               c(0.4000000, 0.2666667, 0.1333333), tolerance = 1e-6)
})

test_that("dtri: density = 0 outside [min, max]", {
  expect_equal(dtri(9,  10, 15, 12), 0)
  expect_equal(dtri(16, 10, 15, 12), 0)
})

test_that("dtri: density >= 0 everywhere", {
  x <- seq(0, 1, by = 0.05)
  expect_true(all(dtri(x) >= 0))
})

test_that("dtri: integrates to 1", {
  # numerical integration over [0,1] with mode=0.5
  x <- seq(0, 1, length.out = 10001)
  dx <- x[2] - x[1]
  integral <- sum(dtri(x)) * dx
  expect_equal(integral, 1, tolerance = 1e-4)
})

test_that("dtri: NA input gives NA output", {
  expect_true(is.na(dtri(NA)))
})

test_that("dtri: mode at x gives peak density 2/(max-min)", {
  # at x=mode density = 2/(max-min)
  expect_equal(dtri(0.5, 0, 1, 0.5), 2, tolerance = 1e-10)
})

# --- ptri ---

test_that("ptri: known values min=2, max=7, mode=5", {
  expect_equal(ptri(3:5, 2, 7, 5),
               c(0.06666667, 0.26666667, 0.60000000), tolerance = 1e-6)
})

test_that("ptri: CDF at min = 0", {
  expect_equal(ptri(0, 0, 1, 0.5), 0)
})

test_that("ptri: CDF at max = 1", {
  expect_equal(ptri(1, 0, 1, 0.5), 1)
})

test_that("ptri: non-decreasing", {
  q <- seq(0, 1, by = 0.1)
  p <- ptri(q)
  expect_true(all(diff(p) >= 0))
})

test_that("ptri: NA input gives NA", {
  expect_true(is.na(ptri(NA)))
})

# --- qTri ---

test_that("qTri: known value min=1, max=4, mode=3 at p=0.25", {
  expect_equal(qTri(0.25, 1, 4, 3), 2.224745, tolerance = 1e-5)
})

test_that("qTri: p=0 returns min", {
  expect_equal(qTri(0, 0, 1, 0.5), 0)
})

test_that("qTri: p=1 returns max", {
  expect_equal(qTri(1, 0, 1, 0.5), 1)
})

test_that("qTri: ptri(qTri(p)) == p roundtrip", {
  p <- c(0.1, 0.25, 0.5, 0.75, 0.9)
  expect_equal(ptri(qTri(p)), p, tolerance = 1e-10)
})

test_that("qTri: invalid params throw error", {
  expect_error(qTri(0.5, 0, 1, 1.5))   # mode > max
  expect_error(qTri(0.5, 0, 1, -0.1))  # mode < min
})


test_that("dtri invalid parameters throw error", {
  expect_error(dtri(0.5, 0, 1, 0))
  expect_error(dtri(0.5, 0, 1, 1))
})

test_that("ptri invalid parameters throw error", {
  expect_error(ptri(0.5, 0, 1, 0))
})

test_that("qTri invalid probability throws error", {
  expect_error(qTri(-0.1))
  expect_error(qTri(1.1))
})

test_that("qTri preserves NA", {
  expect_true(is.na(qTri(NA)))
})

test_that("rTri returns correct length", {
  set.seed(1)
  expect_length(rTri(100), 100)
})

test_that("rTri stays inside support", {
  set.seed(1)
  x <- rTri(1000, 2, 5, 3)
  
  expect_true(all(x >= 2))
  expect_true(all(x <= 5))
})

test_that("rTri invalid n throws error", {
  expect_error(rTri(0))
  expect_error(rTri(-1))
  expect_error(rTri(1.5))
})

