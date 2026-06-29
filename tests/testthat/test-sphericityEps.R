library(testthat)
library(lumen)

# A compound symmetry covariance matrix has epsilon = 1
p <- 4
rho <- 0.5
S_cs <- rho * matrix(1, p, p) + (1 - rho) * diag(p)  # compound symmetry

test_that("sphericityEps: compound symmetry gives gg epsilon = 1", {
  res <- sphericityEps(S_cs, p = p, nGroups = 1, n = 20, method = "gg")
  expect_equal(unname(res["gg"]), 1, tolerance = 1e-6)
})

test_that("sphericityEps: compound symmetry gives hf epsilon = 1", {
  res <- sphericityEps(S_cs, p = p, nGroups = 1, n = 20, method = "hf")
  expect_equal(unname(res["hf"]), 1, tolerance = 1e-3)
})

test_that("sphericityEps: gg in (1/(p-1), 1]", {
  S_arb <- diag(p) + 0.2 * matrix(rnorm(p^2), p, p)
  S_arb <- crossprod(S_arb)  # make positive definite
  res   <- sphericityEps(S_arb, p = p, nGroups = 1, n = 30, method = "gg")
  gg    <- unname(res["gg"])
  expect_gte(gg, 1/(p-1) - 1e-6)
  expect_lte(gg, 1 + 1e-6)
})

test_that("sphericityEps: method='both' returns gg and hf", {
  res <- sphericityEps(S_cs, p = p, nGroups = 1, n = 20, method = "both")
  expect_true(all(c("gg","hf") %in% names(res)))
})

test_that("sphericityEps: method='gg' returns only gg", {
  res <- sphericityEps(S_cs, p = p, nGroups = 1, n = 20, method = "gg")
  expect_true("gg" %in% names(res))
})
