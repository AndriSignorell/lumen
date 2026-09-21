set.seed(1)
x <- rnorm(30, sd = 2)
y <- rnorm(25, sd = 3)

# sample with exactly the requested variance, so that the test statistic
# hits a prescribed value
withVar <- function(n, v, seed = n) {
  set.seed(seed)
  z <- rnorm(n)
  (z - mean(z)) / sd(z) * sqrt(v)
}

# LD reference p-values from an independent implementation
# (scipy.stats, brentq with xtol = rtol = 1e-15)
ldChisq <- data.frame(
  q  = c(15, 3, 12, 1, 30, 4, 0.5, 3),
  df = c(10, 10, 5, 4, 3, 30, 1, 2),
  p  = c(0.16897319272765146, 0.09828410779135308, 0.036137372990705775,
         0.5661292932459128, 1.3800570312945064e-06, 4.002525188555445e-08,
         0.47950012218695337, 0.22313016014842982))

ldF <- data.frame(
  f   = c(2.5, 0.3, 40, 15, 0.1, 3),
  df1 = c(5, 8, 5, 3, 30, 2),
  df2 = c(10, 12, 5, 30, 30, 10),
  p   = c(0.10393193369794569, 0.3835775168416321, 0.0004915823965295844,
          3.8009144117045745e-06, 4.096446692807721e-07, 0.09536743164062497))


# -- htest contract -------------------------------------------------------------

test_that("one-sample: htest contract", {
  r <- varTest(x, sigma2_0 = 4)
  expect_s3_class(r, "htest")
  expect_named(r, c("statistic", "parameter", "p.value", "estimate",
                    "null.value", "alternative", "method", "data.name"))
  expect_identical(names(r$statistic), "X-squared")
  expect_identical(r$parameter, c(df = 29))
  expect_identical(r$estimate, c(variance = var(x)))
  expect_identical(r$null.value, c(variance = 4))
  expect_identical(r$method, "One-sample variance test (classic)")
  expect_identical(r$data.name, "x")
  expect_identical(varTest(x, sigma2_0 = 4, type = "ld")$method,
                   "One-sample variance test (ld)")
})


test_that("two-sample: htest contract", {
  r <- varTest(x, y)
  expect_named(r, c("statistic", "parameter", "p.value", "estimate",
                    "alternative", "method", "data.name"))
  expect_identical(names(r$statistic), "F")
  expect_identical(r$parameter, c(df1 = 29, df2 = 24))
  expect_identical(r$estimate, c("var(x)" = var(x), "var(y)" = var(y)))
  expect_identical(r$method, "Two-sample variance test (classic)")
  expect_identical(r$data.name, "x and y")
})


# -- classic --------------------------------------------------------------------

test_that("classic one-sample equals the chi-squared formulas", {
  q <- 29 * var(x) / 4
  expect_equal(unname(varTest(x, sigma2_0 = 4)$statistic), q)
  expect_equal(varTest(x, sigma2_0 = 4, alternative = "less")$p.value,
               pchisq(q, 29))
  expect_equal(varTest(x, sigma2_0 = 4, alternative = "greater")$p.value,
               pchisq(q, 29, lower.tail = FALSE))
  expect_equal(varTest(x, sigma2_0 = 4)$p.value,
               2 * min(pchisq(q, 29), pchisq(q, 29, lower.tail = FALSE)))
})


test_that("classic two-sample equals var.test() for all alternatives", {
  for (alt in c("two.sided", "less", "greater")) {
    r <- varTest(x, y, alternative = alt)
    v <- var.test(x, y, alternative = alt)
    expect_equal(unname(r$statistic), unname(v$statistic))
    expect_equal(unname(r$parameter), unname(v$parameter))
    expect_equal(r$p.value, v$p.value, label = paste("p-value,", alt))
    expect_identical(r$alternative, alt)
  }
})


test_that("classic: power sanity", {
  set.seed(42)
  expect_gt(varTest(rnorm(200, sd = 2), sigma2_0 = 4)$p.value, 0.05)
  set.seed(42)
  expect_lt(varTest(rnorm(200, sd = 5), sigma2_0 = 1)$p.value, 0.05)
  set.seed(42)
  expect_gt(varTest(rnorm(100, sd = 2), rnorm(100, sd = 2))$p.value, 0.05)
  set.seed(42)
  expect_lt(varTest(rnorm(200, sd = 1), rnorm(200, sd = 5))$p.value, 0.05)
})


# -- lowest-density p-value -----------------------------------------------------

test_that("ld one-sample matches the reference, incl. extreme statistics", {
  for (k in seq_len(nrow(ldChisq))) {
    df <- ldChisq$df[k]
    xk <- withVar(df + 1, ldChisq$q[k] / df)
    r  <- varTest(xk, sigma2_0 = 1, type = "ld")
    expect_equal(unname(r$statistic), ldChisq$q[k], tolerance = 1e-12)
    expect_equal(r$p.value, ldChisq$p[k], tolerance = 1e-8,
                 label = sprintf("LD p, chisq(%g) at %g", df, ldChisq$q[k]))
  }
})


test_that("ld two-sample matches the reference, incl. extreme statistics", {
  for (k in seq_len(nrow(ldF))) {
    xk <- withVar(ldF$df1[k] + 1, ldF$f[k], seed = 1)
    yk <- withVar(ldF$df2[k] + 1, 1, seed = 2)
    r  <- varTest(xk, yk, type = "ld")
    expect_equal(unname(r$statistic), ldF$f[k], tolerance = 1e-12)
    expect_equal(r$p.value, ldF$p[k], tolerance = 1e-8,
                 label = sprintf("LD p, F(%g, %g) at %g",
                                 ldF$df1[k], ldF$df2[k], ldF$f[k]))
  }
})


test_that("ld: one-sided alternatives are the classic ones", {
  for (alt in c("less", "greater")) {
    expect_identical(varTest(x, sigma2_0 = 4, alternative = alt, type = "ld")$p.value,
                     varTest(x, sigma2_0 = 4, alternative = alt)$p.value)
    expect_identical(varTest(x, y, alternative = alt, type = "ld")$p.value,
                     varTest(x, y, alternative = alt)$p.value)
  }
})


test_that(".ldPValue(): both tails have equal density, mode gives p = 1", {
  pd <- function(q, lower.tail = TRUE) pchisq(q, 8, lower.tail = lower.tail)
  dd <- function(q) dchisq(q, 8)
  expect_identical(.ldPValue(6, pd, dd, mode = 6), 1)
  for (q in c(0.5, 3, 11, 25)) {
    p <- .ldPValue(q, pd, dd, mode = 6)
    # the other endpoint: solve P(lower) + P(upper) = p for the matching quantile
    other <- if (q > 6) qchisq(p - pd(q, lower.tail = FALSE), 8) else
      qchisq(p - pd(q), 8, lower.tail = FALSE)
    expect_equal(dd(other), dd(q), tolerance = 1e-7, label = paste("q =", q))
  }
})


test_that("missing values give NA, not an error", {
  expect_identical(varTest(c(x, NA), sigma2_0 = 4)$p.value, NA_real_)
  expect_identical(varTest(c(x, NA), sigma2_0 = 4, type = "ld")$p.value, NA_real_)
  expect_identical(varTest(x, c(y, NA), type = "ld")$p.value, NA_real_)
})


test_that("input validation", {
  expect_error(varTest(x), "sigma2_0 must be provided")
  expect_error(varTest(x, sigma2_0 = 0), "'sigma2_0' must be")
  expect_error(varTest(x, sigma2_0 = -1), "'sigma2_0' must be")
  expect_error(varTest(x, sigma2_0 = c(1, 2)), "'sigma2_0' must be")
  expect_error(varTest(1, sigma2_0 = 1), "at least two")
  expect_error(varTest(letters, sigma2_0 = 1), "'x' must be")
  expect_error(varTest(x, 1), "'y' must be")
  expect_error(varTest(x, sigma2_0 = 4, alternative = "both"))
  expect_error(varTest(x, sigma2_0 = 4, type = "exact"))
})


# -- formula interface ----------------------------------------------------------

test_that("formula method equals the default method", {
  d <- data.frame(v = c(x, y), g = factor(rep(c("A", "B"), c(30, 25))))
  r <- varTest(v ~ g, data = d, type = "ld", alternative = "greater")
  s <- varTest(x, y, type = "ld", alternative = "greater")
  expect_s3_class(r, "htest")
  expect_equal(r[c("statistic", "parameter", "p.value", "estimate")],
               s[c("statistic", "parameter", "p.value", "estimate")])
  expect_type(r$data.name, "character")
  expect_length(r$data.name, 1L)
})


test_that("formula method: subset is evaluated in data", {
  d <- data.frame(v = c(x, y), g = factor(rep(c("A", "B"), c(30, 25))),
                  keep = rep(c(TRUE, FALSE, TRUE), c(20, 10, 25)))
  r <- varTest(v ~ g, data = d, subset = keep)
  s <- varTest(x[1:20], y)
  expect_equal(r$statistic, s$statistic)
  expect_identical(r$parameter, c(df1 = 19, df2 = 24))
})


test_that("formula method: errors", {
  d <- data.frame(v = 1:9, g = factor(rep(c("A", "B", "C"), 3)))
  expect_error(varTest(v ~ g, data = d))
  expect_error(varTest.formula(~ g, data = d), "'formula' missing or incorrect")
})
