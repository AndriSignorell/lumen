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

# --- qtri ---

test_that("qtri: known value min=1, max=4, mode=3 at p=0.25", {
  expect_equal(qtri(0.25, 1, 4, 3), 2.224745, tolerance = 1e-5)
})

test_that("qtri: p=0 returns min", {
  expect_equal(qtri(0, 0, 1, 0.5), 0)
})

test_that("qtri: p=1 returns max", {
  expect_equal(qtri(1, 0, 1, 0.5), 1)
})

test_that("qtri: ptri(qtri(p)) == p roundtrip", {
  p <- c(0.1, 0.25, 0.5, 0.75, 0.9)
  expect_equal(ptri(qtri(p)), p, tolerance = 1e-10)
})

test_that("qtri: invalid params throw error", {
  expect_error(qtri(0.5, 0, 1, 1.5))   # mode > max
  expect_error(qtri(0.5, 0, 1, -0.1))  # mode < min
})


test_that("dtri invalid parameters throw error", {
  expect_error(dtri(0.5, 0, 1, 0))
  expect_error(dtri(0.5, 0, 1, 1))
})

test_that("ptri invalid parameters throw error", {
  expect_error(ptri(0.5, 0, 1, 0))
})

test_that("qtri: p outside [0,1] gives NaN with a warning", {
  expect_warning(res <- qtri(c(-0.1, 1.1)), "NaN")
  expect_true(all(is.nan(res)))
})

test_that("dtri, ptri and qtri handle log, lower.tail and log.p", {
  x <- c(2.5, 3, 5, 6)
  expect_equal(dtri(x, 2, 7, 5, log = TRUE), log(dtri(x, 2, 7, 5)))
  expect_equal(ptri(x, 2, 7, 5, log.p = TRUE), log(ptri(x, 2, 7, 5)))
  expect_equal(ptri(x, 2, 7, 5, lower.tail = FALSE), 1 - ptri(x, 2, 7, 5))
  p <- c(0.1, 0.5, 0.9)
  expect_equal(qtri(log(p), 2, 7, 5, log.p = TRUE), qtri(p, 2, 7, 5))
  expect_equal(qtri(p, 2, 7, 5, lower.tail = FALSE), qtri(1 - p, 2, 7, 5))
})

test_that("dtri and ptri recycle their parameters and keep names", {
  expect_equal(dtri(c(0.25, 0.5), min = 0, max = c(1, 2), mode = c(0.5, 1)),
               c(dtri(0.25, 0, 1, 0.5), dtri(0.5, 0, 2, 1)))
  expect_named(dtri(c(a = 0.25, b = 0.5)), c("a", "b"))
  expect_named(ptri(c(a = 0.25, b = 0.5)), c("a", "b"))
})

test_that("mtri agrees with a numerical moment of dtri", {
  m <- mtri(2, 7, 5)
  mu <- integrate(function(x) x * dtri(x, 2, 7, 5), 2, 7)$value
  v  <- integrate(function(x) (x - mu)^2 * dtri(x, 2, 7, 5), 2, 7)$value
  expect_equal(unname(m["mean"]),     mu, tolerance = 1e-4)
  expect_equal(unname(m["variance"]), v,  tolerance = 1e-4)
})

test_that("qtri preserves NA", {
  expect_true(is.na(qtri(NA)))
})

test_that("rtri returns correct length", {
  set.seed(1)
  expect_length(rtri(100), 100)
})

test_that("rtri stays inside support", {
  set.seed(1)
  x <- rtri(1000, 2, 5, 3)
  
  expect_true(all(x >= 2))
  expect_true(all(x <= 5))
})

test_that("rtri invalid n throws error", {
  expect_error(rtri(0))
  expect_error(rtri(-1))
  expect_error(rtri(1.5))
})

