# Blaker interval: blaker_find_crossing_cpp() and .binomCI.blaker() ----------

# acceptability straight from Blaker's definition, P(gamma(X) <= gamma(x)),
# gamma(k) = min(P(X >= k), P(X <= k)); independent of the C++ tail algebra
blakerAccept <- function(x, n, p) {
  f  <- dbinom(0:n, n, p)
  lo <- cumsum(f)
  up <- rev(cumsum(rev(f)))
  g  <- pmin(lo, up)
  sum(f[g <= g[x + 1L] * (1 + 1e-7)])
}

bfc <- function(...) lumen:::blaker_find_crossing_cpp(...)

test_that(".binomCI.blaker limits sit on the acceptability crossing", {
  alpha <- 0.05
  for (xn in list(c(7, 20), c(1, 10), c(9, 10), c(3, 10000),
                  c(50000, 100000))) {
    x <- xn[1]; n <- xn[2]
    ci <- lumen:::.binomCI.blaker(x, n, alpha)
    lci <- ci[["lci"]]; uci <- ci[["uci"]]
    cp <- c(qbeta(alpha / 2, x, n - x + 1), qbeta(1 - alpha / 2, x + 1, n - x))
    # well above the bisection tolerance (1e-12 relative), well below the
    # distance to the next jump of the acceptability function
    el <- 1e-6 * lci
    eu <- 1e-6 * uci

    info <- sprintf("x = %d, n = %d", x, n)
    expect_true(lci < x / n && x / n < uci, info = info)
    # inside Clopper-Pearson
    expect_true(lci >= cp[1] && uci <= cp[2], info = info)
    # crossing: acceptable at the limit, not acceptable just outside
    expect_gte(blakerAccept(x, n, lci + el), alpha - 1e-9)
    expect_lt(blakerAccept(x, n, lci - el), alpha)
    expect_gte(blakerAccept(x, n, uci - eu), alpha - 1e-9)
    expect_lt(blakerAccept(x, n, uci + eu), alpha)
  }
})

test_that(".binomCI.blaker boundary counts", {
  expect_equal(lumen:::.binomCI.blaker(0, 10, 0.05)[["lci"]], 0)
  expect_equal(lumen:::.binomCI.blaker(10, 10, 0.05)[["uci"]], 1)
  u <- lumen:::.binomCI.blaker(0, 10, 0.05)[["uci"]]
  expect_gte(blakerAccept(0, 10, u * (1 - 1e-7)), 0.05 - 1e-9)
  expect_lt(blakerAccept(0, 10, u * (1 + 1e-7)), 0.05)
})

test_that("blaker_find_crossing_cpp: degenerate bracket", {
  expect_identical(bfc(5, 10, 0.05, 0.5, 0.5, TRUE), 0.5)
  expect_identical(bfc(5, 10, 0.05, 0.6, 0.4, TRUE), 0.4)
  expect_identical(bfc(5, 10, 0.05, 0.6, 0.4, FALSE), 0.6)
})

test_that("blaker_find_crossing_cpp: endpoint already acceptable", {
  # accept_bin(x, n, x/n) = 1
  expect_identical(bfc(5, 10, 0.05, 0.5, 0.9, TRUE), 0.5)
  expect_identical(bfc(5, 10, 0.05, 0.1, 0.5, FALSE), 0.5)
  # p = 0 resp. p = 1 are acceptable only for x = 0 resp. x = n
  expect_identical(bfc(0, 10, 0.05, 0, 0.5, TRUE), 0)
  expect_identical(bfc(10, 10, 0.05, 0.5, 1, FALSE), 1)
})

test_that("blaker_find_crossing_cpp: coarse scan recovers an unbracketed crossing", {
  ref <- lumen:::.binomCI.blaker(5, 10, 0.05)
  # [0, 1]: both ends below alpha, the scan has to locate the window
  expect_equal(bfc(5, 10, 0.05, 0, 1, TRUE),  ref[["lci"]], tolerance = 1e-10)
  expect_equal(bfc(5, 10, 0.05, 0, 1, FALSE), ref[["uci"]], tolerance = 1e-10)
})

test_that("blaker_find_crossing_cpp: conservative fallback when nothing is found", {
  # bracket entirely outside the acceptance region
  expect_identical(bfc(5, 10, 0.05, 0, 1e-6, TRUE), 1e-6)
  expect_identical(bfc(5, 10, 0.05, 1 - 1e-6, 1, FALSE), 1 - 1e-6)
})

test_that("blaker_find_crossing_cpp: tol = 0 stops at machine precision", {
  ref <- lumen:::.binomCI.blaker(7, 20, 0.05)
  cp <- qbeta(0.025, 7, 14)
  expect_equal(bfc(7, 20, 0.05, cp, 7 / 20, TRUE, tol = 0), ref[["lci"]],
               tolerance = 1e-11)
})
