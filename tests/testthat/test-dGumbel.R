library(testthat)
library(lumen)

tol <- 1e-6

# --- dgumbel ---

test_that("dgumbel: density >= 0", {
  expect_true(all(dgumbel(seq(-5, 5, by = 0.5)) >= 0))
})

test_that("dgumbel: integrates to 1", {
  x <- seq(-10, 30, length.out = 100001)
  dx <- x[2] - x[1]
  expect_equal(sum(dgumbel(x, loc = 2, scale = 3)) * dx, 1, tolerance = 1e-4)
})

test_that("dgumbel: log=TRUE returns log of density", {
  x <- c(-1, 0, 1, 2)
  expect_equal(dgumbel(x, log = TRUE), log(dgumbel(x)), tolerance = tol)
})

test_that("dgumbel: mode at loc (density maximum)", {
  # mode of Gumbel = loc; density at loc = 1/(scale*e)
  expect_equal(dgumbel(0, loc = 0, scale = 1), exp(-1), tolerance = tol)
})

# --- pgumbel ---

test_that("pgumbel: in [0,1]", {
  expect_true(all(pgumbel(seq(-5, 10)) >= 0))
  expect_true(all(pgumbel(seq(-5, 10)) <= 1))
})

test_that("pgumbel: non-decreasing", {
  q <- seq(-5, 10, by = 0.5)
  expect_true(all(diff(pgumbel(q)) >= 0))
})

test_that("pgumbel: lower.tail=FALSE complement", {
  q <- c(-2, 0, 2)
  expect_equal(pgumbel(q, lower.tail = FALSE), 1 - pgumbel(q), tolerance = tol)
})

# --- qgumbel ---

test_that("qgumbel: pgumbel(qgumbel(p)) roundtrip", {
  p <- c(0.1, 0.25, 0.5, 0.75, 0.9)
  expect_equal(pgumbel(qgumbel(p, loc = 1, scale = 2), loc = 1, scale = 2),
               p, tolerance = tol)
})

test_that("qgumbel: median = loc - scale*log(log(2))", {
  # F^{-1}(0.5) = loc - scale * log(-log(0.5))
  expect_equal(qgumbel(0.5, loc = 2, scale = 1),
               2 - log(-log(0.5)), tolerance = tol)
})

test_that("qgumbel: invalid p throws error", {
  expect_error(qgumbel(0))
  expect_error(qgumbel(1))
})

# --- rgumbel ---

test_that("rgumbel: returns correct length", {
  set.seed(1)
  expect_length(rgumbel(100), 100)
})

test_that("rgumbel: sample mean ~ loc + scale * Euler_gamma", {
  set.seed(42)
  r <- rgumbel(50000, loc = 2, scale = 1)
  euler <- 0.5772156649
  expect_equal(mean(r), 2 + euler, tolerance = 0.02)
})


test_that("pgumbel lower.tail roundtrip", {
  q <- c(-2, 0, 2)
  
  expect_equal(
    pgumbel(q, lower.tail = FALSE),
    1 - pgumbel(q)
  )
})

test_that("qgumbel lower.tail=FALSE roundtrip", {
  p <- c(0.1, 0.5, 0.9)
  
  expect_equal(
    pgumbel(
      qgumbel(p, lower.tail = FALSE),
      lower.tail = FALSE
    ),
    p,
    tolerance = 1e-6
  )
})

test_that("dgumbel invalid scale throws error", {
  expect_error(dgumbel(1, scale = -1))
})

test_that("rgumbel invalid scale throws error", {
  expect_error(rgumbel(10, scale = -1))
})

