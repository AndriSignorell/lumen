# Null distribution of the Anderson-Darling statistic: the Rcpp kernels in
# src/ad_test.cpp and the R wrappers pAD()/qAD().
#
# goftest ships the same Marsaglia C code and is used as reference, EXCEPT
# on the z-intervals where ADf() evaluates cPhi() at |x| >= 17
# (z in [(4j+1)^2 * 1.2337 / 150, (4j+1)^2 * 1.2337 / 144.5], j = 1, 2, ...):
# there the original code reads past its lookup table (goftest 1.2-3 returns
# NaN for pAD(0.21, fast = FALSE)). Those intervals are covered by the
# exact-vs-approximation test below instead.

qSafe <- c(0.05, 0.3, 0.5, 1, 1.99, 2, 2.01, 3, 4.5, 8, 15)   # outside them
qHit  <- c(0.21, 0.68, 1.42, 2.42, 3.7, 5.2, 7.0, 9.1)        # inside them


test_that("C++ kernels agree with goftest", {
  skip_if_not_installed("goftest")

  expect_equal(ad_prob_exact_inf_cpp(qSafe),
               goftest::pAD(qSafe, fast = FALSE), tolerance = 1e-12)
  q <- c(qSafe, qHit)
  expect_equal(ad_prob_approx_inf_cpp(q),
               goftest::pAD(q, fast = TRUE), tolerance = 1e-12)
  # finite n: q >= 0.4 because goftest returns values < 0 below that
  q <- q[q >= 0.4]
  for (n in c(1L, 5L, 20L, 200L))
    expect_equal(ad_prob_n_cpp(q, n), goftest::pAD(q, n = n),
                 tolerance = 1e-12, label = paste("n =", n))
})


test_that("exact and approximate asymptotic cdf agree on a dense grid", {
  # Marsaglia's adinf() is accurate to ~2e-5; a garbage read in cPhi() on
  # the qHit intervals would show up here
  z <- sort(c(seq(0.02, 30, by = 0.005), qHit))
  pe <- ad_prob_exact_inf_cpp(z)
  expect_true(all(is.finite(pe)))
  expect_lt(max(abs(pe - ad_prob_approx_inf_cpp(z))), 5e-5)
  expect_true(all(diff(pe) > -1e-12))          # monotone up to rounding
})


test_that("asymptotic cdf reproduces the classical critical values", {
  # upper 10% and 5% points of A^2, all parameters known (Stephens 1974)
  expect_equal(ad_prob_exact_inf_cpp(c(1.933, 2.492)), c(0.90, 0.95),
               tolerance = 1e-4)
  expect_equal(ad_prob_approx_inf_cpp(c(1.933, 2.492)), c(0.90, 0.95),
               tolerance = 1e-4)
})


test_that("finite-n cdf converges to the asymptotic one and stays in [0, 1]", {
  z <- seq(0.001, 40, by = 0.01)
  expect_lt(max(abs(ad_prob_n_cpp(z, 10000L) - ad_prob_exact_inf_cpp(z))), 5e-5)
  for (n in c(1L, 2L, 5L, 10L, 50L)) {
    p <- ad_prob_n_cpp(z, n)
    expect_true(all(p >= 0 & p <= 1), label = paste("range for n =", n))
  }
})


test_that("C++ kernels handle boundary and non-finite input", {
  b <- c(-1, 0, Inf)
  expect_identical(ad_prob_exact_inf_cpp(b), c(0, 0, 1))
  expect_identical(ad_prob_approx_inf_cpp(b), c(0, 0, 1))
  expect_identical(ad_prob_n_cpp(b, 5L), c(0, 0, 1))
  # NaN used to crash the session in cPhi() (int cast of NaN as array index)
  expect_true(is.na(ad_prob_exact_inf_cpp(NaN)))
  expect_true(is.na(ad_prob_approx_inf_cpp(NaN)))
  expect_true(is.na(ad_prob_n_cpp(NaN, 5L)))
  expect_length(ad_prob_exact_inf_cpp(numeric(0)), 0L)
})


# -- R wrappers (goftest-compatible signature) --------------------------------

test_that("pAD() matches goftest::pAD()", {
  skip_if_not_installed("goftest")
  # finite n only from q = 0.4 up: below that goftest returns values outside
  # [0, 1] (pAD(0.05, n = 10) = -1.3e-7), see the clamping test below
  for (n in c(Inf, 10)) {
    q <- c(-1, 0, if (is.finite(n)) qSafe[qSafe >= 0.4] else qSafe, Inf)
    for (fast in c(TRUE, FALSE)) for (lt in c(TRUE, FALSE))
      expect_equal(pAD(q, n = n, lower.tail = lt, fast = fast),
                   goftest::pAD(q, n = n, lower.tail = lt, fast = fast),
                   tolerance = 1e-12,
                   label = sprintf("n = %s, fast = %s, lower.tail = %s", n, fast, lt))
  }
})


test_that("pAD() clamps where Marsaglia's correction leaves [0, 1]", {
  # errfix() pushes the cdf below 0 for small q at small n -- goftest 1.2-3
  # passes that through, lumen does not
  q <- seq(0.001, 0.5, by = 0.001)
  for (n in c(1L, 2L, 5L, 10L, 50L)) {
    p <- pAD(q, n = n)
    expect_true(all(p >= 0 & p <= 1), label = paste("range for n =", n))
    expect_true(all(pAD(q, n = n, lower.tail = FALSE) <= 1),
                label = paste("upper tail for n =", n))
  }
})


test_that("pAD(): lower.tail, boundaries, vectorisation", {
  q <- c(0.5, 1, 2.5)
  expect_equal(pAD(q, lower.tail = FALSE), 1 - pAD(q))
  expect_equal(pAD(q, n = 12, lower.tail = FALSE), 1 - pAD(q, n = 12))
  expect_identical(pAD(c(-Inf, -1, 0, Inf)), c(0, 0, 0, 1))
  expect_length(pAD(numeric(0)), 0L)
})


test_that("qAD() inverts pAD()", {
  p <- c(0.1, 0.5, 0.9, 0.95, 0.99)
  expect_equal(pAD(qAD(p)), p, tolerance = 1e-4)
  expect_equal(pAD(qAD(p, n = 20), n = 20), p, tolerance = 1e-4)
  expect_equal(qAD(0.95), 2.492, tolerance = 1e-3)
})
