library(testthat)
library(lumen)

tol <- 1e-6

test_that("dRevGumbel: density >= 0", {
  expect_true(all(dRevGumbel(seq(-3, 3, by = 0.25)) >= 0))
})

test_that("dRevGumbel: integrates to 1", {
  x <- seq(-20, 10, length.out = 100001)
  dx <- x[2] - x[1]
  expect_equal(sum(dRevGumbel(x, location = 0, scale = 1)) * dx, 1,
               tolerance = 1e-4)
})

test_that("dRevGumbel: mode at location", {
  # mode of RevGumbel = location; d/dx = 0 at x = location + scale * log(1) = location
  # density at mode = exp(0)*exp(-1)/scale = 1/(e*scale)
  expect_equal(dRevGumbel(0, location = 0, scale = 1), exp(-1), tolerance = tol)
})

test_that("dRevGumbel: invalid scale throws error", {
  expect_error(dRevGumbel(1, scale = -1))
  expect_error(dRevGumbel(1, scale = 0))
})

test_that("pRevGumbel: in [0,1]", {
  q <- seq(-5, 5, by = 0.5)
  p <- pRevGumbel(q)
  expect_true(all(p >= 0 & p <= 1))
})

test_that("pRevGumbel: non-decreasing", {
  q <- seq(-5, 5, by = 0.25)
  expect_true(all(diff(pRevGumbel(q)) >= 0))
})

test_that("pRevGumbel: relation to pgumbel", {
  # RevGumbel(loc, scale) = -Gumbel(-loc, scale)
  # pRevGumbel(q) = 1 - pgumbel(-q, loc=-loc, scale=scale)
  q <- c(-2, -1, 0, 1)
  expect_equal(pRevGumbel(q, location = 0, scale = 1),
               1 - pgumbel(-q, loc = 0, scale = 1), tolerance = tol)
})

test_that("qRevGumbel: pRevGumbel(qRevGumbel(1-p)) == p roundtrip", {
  # qRevGumbel(p) = loc + scale*log(-log(p)) inverts the survival function:
  # 1 - pRevGumbel(q) = exp(-exp((q-loc)/scale))
  # => pRevGumbel(qRevGumbel(1-p)) == p
  p <- c(0.1, 0.25, 0.5, 0.75, 0.9)
  expect_equal(pRevGumbel(qRevGumbel(1 - p, location = 1, scale = 2),
                          location = 1, scale = 2), p, tolerance = tol)
})

test_that("rRevGumbel: returns correct length", {
  set.seed(1)
  expect_length(rRevGumbel(50), 50)
})


test_that("pRevGumbel invalid scale throws error", {
  expect_error(pRevGumbel(1, scale = 0))
})

test_that("qRevGumbel invalid scale throws error", {
  expect_error(qRevGumbel(0.5, scale = 0))
})

test_that("qRevGumbelExp equals exp(qRevGumbel())", {
  p <- c(0.2, 0.5, 0.8)
  
  expect_equal(
    qRevGumbelExp(p),
    exp(qRevGumbel(p))
  )
})

test_that("rRevGumbel invalid scale throws error", {
  expect_error(rRevGumbel(10, scale = -1))
})

