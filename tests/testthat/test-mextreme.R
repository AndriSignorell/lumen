library(testthat)
library(lumen)

tol <- 1e-8

# The moments are checked against numerical integration of the corresponding
# density, which ties the m*() functions to the d*() functions they describe.

# integrate() probes the infinite ends of the range, where x * d(x) would be
# Inf * 0; the integrand is zero there and is set so explicitly
mNumeric <- function(dFun, lower = -Inf, upper = Inf, ...) {
  wrap <- function(g) function(x) ifelse(is.finite(x), g(x), 0)
  mu <- integrate(wrap(function(x) x * dFun(x, ...)), lower, upper)$value
  v  <- integrate(wrap(function(x) (x - mu)^2 * dFun(x, ...)),
                  lower, upper)$value
  c(mean = mu, variance = v)
}

# --- mgumbel ---

test_that("mgumbel: mean = loc + scale * gamma, variance = pi^2/6 * scale^2", {
  r <- mgumbel(loc = 2, scale = 3)
  expect_equal(unname(r["mean"]),     2 + 3 * -digamma(1), tolerance = tol)
  expect_equal(unname(r["variance"]), pi^2 / 6 * 9,        tolerance = tol)
})

test_that("mgumbel agrees with the numerical moments of dgumbel", {
  expect_equal(mgumbel(2, 3), mNumeric(dgumbel, loc = 2, scale = 3),
               tolerance = 1e-6)
})

# --- mrevgumbel ---

test_that("mrevgumbel: the reflection flips the sign of the scale term", {
  r <- mrevgumbel(loc = 2, scale = 3)
  expect_equal(unname(r["mean"]),     2 - 3 * -digamma(1), tolerance = tol)
  expect_equal(unname(r["variance"]), pi^2 / 6 * 9,        tolerance = tol)
})

test_that("mrevgumbel agrees with the numerical moments of drevgumbel", {
  expect_equal(mrevgumbel(2, 3),
               mNumeric(drevgumbel, loc = 2, scale = 3),
               tolerance = 1e-6)
})

test_that("mrevgumbel mirrors mgumbel about loc", {
  expect_equal(unname(mrevgumbel(0, 2)["mean"]), -unname(mgumbel(0, 2)["mean"]))
  expect_equal(mrevgumbel(0, 2)["variance"], mgumbel(0, 2)["variance"])
})

# --- mfrechet ---

test_that("mfrechet: moments exist only for shape > 1 resp. > 2", {
  expect_true(is.na(mfrechet(shape = 0.5)["mean"]))
  expect_true(is.na(mfrechet(shape = 1.5)["variance"]))
  expect_false(is.na(mfrechet(shape = 2.5)["variance"]))
})

test_that("mfrechet agrees with the numerical moments of dfrechet", {
  expect_equal(mfrechet(0, 1, 3), mNumeric(dfrechet, 0, Inf, loc = 0,
                                           scale = 1, shape = 3),
               tolerance = 1e-4)
})

# --- mrevweibull ---

test_that("mrevweibull agrees with the numerical moments of drevweibull", {
  expect_equal(mrevweibull(0, 1, 2), mNumeric(drevweibull, -Inf, 0, loc = 0,
                                              scale = 1, shape = 2),
               tolerance = 1e-6)
})

# --- mgev ---

test_that("mgev: shape = 0 reproduces the Gumbel moments", {
  expect_equal(mgev(2, 3, 0), mgumbel(2, 3), tolerance = tol)
})

test_that("mgev: moments exist only for shape < 1 resp. < 1/2", {
  expect_true(is.na(mgev(shape = 1.5)["mean"]))
  expect_true(is.na(mgev(shape = 0.7)["variance"]))
  expect_false(is.na(mgev(shape = 0.3)["variance"]))
})

test_that("mgev agrees with the numerical moments of dgev", {
  # shape > 0 is bounded below at loc - scale/shape, shape < 0 above
  expect_equal(mgev(1, 2, 0.3), mNumeric(dgev, -5.7, Inf, loc = 1, scale = 2,
                                         shape = 0.3),
               tolerance = 1e-4)
  expect_equal(mgev(1, 2, -0.3), mNumeric(dgev, -Inf, 7.7, loc = 1, scale = 2,
                                          shape = -0.3),
               tolerance = 1e-5)
})

# --- mgpd ---

test_that("mgpd: moments exist only for shape < 1 resp. < 1/2", {
  expect_true(is.na(mgpd(shape = 1.5)["mean"]))
  expect_true(is.na(mgpd(shape = 0.7)["variance"]))
})

test_that("mgpd: shape = 0 is the exponential shifted by loc", {
  expect_equal(unname(mgpd(loc = 1, scale = 2, shape = 0)),
               unname(c(1 + 2, 4)), tolerance = tol)
})

test_that("mgpd agrees with the numerical moments of dgpd", {
  expect_equal(mgpd(1, 2, 0.3), mNumeric(dgpd, 1, Inf, loc = 1, scale = 2,
                                         shape = 0.3),
               tolerance = 1e-4)
})

# --- mgompertz ---

test_that("mgompertz: shape = 0 gives the exponential moments", {
  expect_equal(unname(mgompertz(shape = 0, rate = 2)), c(0.5, 0.25),
               tolerance = tol)
})

test_that("mgompertz: a negative shape leaves the moments undefined", {
  expect_true(all(is.na(mgompertz(shape = -0.5, rate = 1))))
})

test_that("mgompertz agrees with the numerical moments of dgompertz", {
  expect_equal(mgompertz(shape = 1, rate = 1),
               # a finite upper limit: the density is nil well before
               # exp(shape * x) overflows
               mNumeric(dgompertz, 0, 100, shape = 1, rate = 1),
               tolerance = 1e-6)
})

# --- argument validation, shared by all of them ---

test_that("the moment functions reject invalid parameters", {
  expect_error(mgumbel(scale = 0))
  expect_error(mrevgumbel(scale = -1))
  expect_error(mfrechet(shape = 0))
  expect_error(mrevweibull(scale = 0))
  expect_error(mgev(scale = -1))
  expect_error(mgpd(shape = c(0, 1)))
  expect_error(mgompertz(shape = 1, rate = 0))
  expect_error(mgumbel(loc = NA))
})
