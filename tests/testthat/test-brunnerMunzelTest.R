
# Reference values were produced from the Brunner & Munzel (2000) pain-score
# data, the example distributed with lawstat::brunner.munzel.test. Exact
# permutation counts were cross-checked against an independent brute-force
# enumeration.

x <- c(1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 2, 4, 1, 1)
y <- c(3, 3, 4, 3, 1, 2, 3, 1, 1, 5, 4)


test_that("brunnerMunzelTest reproduces the published reference values", {

  r <- brunnerMunzelTest(x, y)

  expect_equal(unname(r$statistic), 3.13746748230295, tolerance = 1e-12)
  expect_equal(unname(r$parameter), 17.6828419794815, tolerance = 1e-12)
  expect_equal(unname(r$estimate),  0.788961038961039, tolerance = 1e-12)
  expect_equal(r$stderr,            0.0921000904681686, tolerance = 1e-12)
  expect_equal(r$p.value,           0.00578620866615148, tolerance = 1e-12)
  expect_equal(as.vector(r$conf.int),
               c(0.595216864253736, 0.982705213668342), tolerance = 1e-12)

  expect_s3_class(r, "htest")
  expect_identical(names(r$statistic), "W")
  expect_identical(names(r$parameter), "df")
  expect_identical(names(r$estimate), "P(X < Y) + 0.5 * P(X = Y)")
  expect_equal(attr(r$conf.int, "conf.level"), 0.95)
})


test_that("the estimate equals the Mann-Whitney common language effect size", {

  set.seed(11)
  a <- rnorm(13)
  b <- rnorm(17, mean = 0.6, sd = 2.5)

  u <- mean(outer(a, b, "<")) + 0.5 * mean(outer(a, b, "=="))
  expect_equal(unname(brunnerMunzelTest(a, b)$estimate), u)

  # ... and with ties
  a <- c(1, 1, 1, 2, 2, 3, 3, 3, 3, 4)
  b <- c(2, 2, 3, 3, 3, 4, 4, 5, 5, 5)
  u <- mean(outer(a, b, "<")) + 0.5 * mean(outer(a, b, "=="))
  expect_equal(unname(brunnerMunzelTest(a, b)$estimate), u)
  expect_equal(unname(brunnerMunzelTest(a, b)$statistic), 2.80467836583555,
               tolerance = 1e-12)
  expect_equal(unname(brunnerMunzelTest(a, b)$parameter), 17.6571100442069,
               tolerance = 1e-12)
})


test_that("swapping the samples mirrors the result", {

  f <- brunnerMunzelTest(x, y)
  g <- brunnerMunzelTest(y, x)

  expect_equal(unname(f$estimate) + unname(g$estimate), 1)
  expect_equal(unname(f$statistic), -unname(g$statistic))
  expect_equal(unname(f$parameter), unname(g$parameter))
  expect_equal(f$p.value, g$p.value)
  expect_equal(f$stderr, g$stderr)

  # one-sided alternatives swap accordingly
  expect_equal(brunnerMunzelTest(x, y, alternative = "greater")$p.value,
               brunnerMunzelTest(y, x, alternative = "less")$p.value)
})


test_that("the test is invariant under monotone transformation of the data", {

  f <- brunnerMunzelTest(x, y)
  g <- brunnerMunzelTest(exp(x), exp(y))
  expect_equal(unname(f$statistic), unname(g$statistic))
  expect_equal(unname(f$estimate), unname(g$estimate))
})


test_that("alternative is stated in terms of p, not of x versus y", {

  # y is stochastically larger, so "greater" must be the significant direction
  expect_equal(brunnerMunzelTest(x, y, alternative = "greater")$p.value,
               0.00289310433307574, tolerance = 1e-12)
  expect_equal(brunnerMunzelTest(x, y, alternative = "less")$p.value,
               0.997106895666924, tolerance = 1e-12)

  # one-sided intervals report the open end at the range limit, not +/-Inf
  ci <- brunnerMunzelTest(x, y, alternative = "greater")$conf.int
  expect_equal(as.vector(ci), c(0.629098293325597, 1), tolerance = 1e-12)
  ci <- brunnerMunzelTest(x, y, alternative = "less")$conf.int
  expect_equal(as.vector(ci), c(0, 0.948823784596481), tolerance = 1e-12)
})


test_that("method = 'normal' drops the Satterthwaite correction", {

  r <- brunnerMunzelTest(x, y, method = "normal")
  expect_equal(unname(r$parameter), Inf)
  expect_equal(r$p.value, 0.00170414176003821, tolerance = 1e-12)
  expect_equal(r$p.value, 2 * pnorm(-abs(unname(r$statistic))))
})


test_that("p0 shifts the null value", {

  r <- brunnerMunzelTest(x, y, p0 = 0.7)
  expect_equal(unname(r$statistic), 0.965916955225853, tolerance = 1e-12)
  expect_equal(r$p.value, 0.347111411309979, tolerance = 1e-12)
  expect_equal(unname(r$null.value), 0.7)
  # the estimate and its interval do not depend on p0
  expect_equal(unname(r$estimate), unname(brunnerMunzelTest(x, y)$estimate))
  expect_equal(r$conf.int, brunnerMunzelTest(x, y)$conf.int)
})


test_that("the exact permutation test matches brute-force enumeration", {

  # C(12, 5) = 792 splits, counts verified against a full independent enumeration
  a <- c(1, 1, 2, 3, 3)
  b <- c(4, 2, 5, 1, 6, 3, 7)

  expect_equal(brunnerMunzelTest(a, b, method = "permutation")$p.value,
               75 / 792, tolerance = 1e-12)
  expect_equal(brunnerMunzelTest(a, b, method = "permutation",
                                 alternative = "greater")$p.value,
               44 / 792, tolerance = 1e-12)
  expect_equal(brunnerMunzelTest(a, b, method = "permutation",
                                 alternative = "less")$p.value,
               766 / 792, tolerance = 1e-12)

  r <- brunnerMunzelTest(a, b, method = "permutation")
  expect_match(r$method, "exact studentized permutation, 792 splits", fixed = TRUE)
  expect_null(r$parameter)
  # the reported statistic is still the studentized one
  expect_equal(unname(r$statistic), 2.09656967344384, tolerance = 1e-12)
})


test_that("the exact permutation p-value is invariant to input order", {

  a <- c(1, 1, 2, 3, 3)
  b <- c(4, 2, 5, 1, 6, 3, 7)
  set.seed(3)
  expect_equal(brunnerMunzelTest(sample(a), sample(b), method = "permutation")$p.value,
               brunnerMunzelTest(a, b, method = "permutation")$p.value)
})


test_that("exact and Monte-Carlo permutation agree to sampling error", {

  a <- c(1, 1, 2, 3, 3)
  b <- c(4, 2, 5, 1, 6, 3, 7)

  set.seed(42)
  mc <- brunnerMunzelTest(a, b, method = "permutation", exact = FALSE,
                          nPerm = 2e5)$p.value
  expect_equal(mc, 75 / 792, tolerance = 0.01)
  expect_match(brunnerMunzelTest(a, b, method = "permutation", exact = FALSE,
                                 nPerm = 100)$method, "100 resamples", fixed = TRUE)
})


test_that("Monte-Carlo p-values are never zero and are reproducible", {

  set.seed(1)
  p1 <- brunnerMunzelTest(x, y, method = "permutation", exact = FALSE,
                          nPerm = 500)$p.value
  set.seed(1)
  p2 <- brunnerMunzelTest(x, y, method = "permutation", exact = FALSE,
                          nPerm = 500)$p.value

  expect_identical(p1, p2)
  expect_gte(p1, 1 / 501)
})


test_that("exact enumeration is used only below the split threshold", {

  # C(25, 14) = 4'457'400 exceeds the default cutoff of 1e6
  set.seed(5)
  expect_match(brunnerMunzelTest(x, y, method = "permutation")$method,
               "resamples", fixed = TRUE)
  expect_match(brunnerMunzelTest(x, y, method = "permutation", exact = TRUE)$method,
               "4457400 splits", fixed = TRUE)
  expect_equal(brunnerMunzelTest(x, y, method = "permutation", exact = TRUE)$p.value,
               0.00803764526405528, tolerance = 1e-12)
})


test_that("non-overlapping samples fall back to the exact permutation test", {

  expect_warning(r <- brunnerMunzelTest(1:3, 4:6), "variance estimator is zero")

  # smallest attainable two-sided p-value is 2 / C(6, 3) = 0.1
  expect_equal(r$p.value, 0.1)
  expect_equal(unname(r$estimate), 1)
  expect_equal(unname(r$statistic), Inf)
  expect_match(r$method, "exact studentized permutation", fixed = TRUE)

  # a zero Wald standard error is an empty variance estimate, not certainty
  expect_equal(as.vector(r$conf.int), c(NA_real_, NA_real_))
  expect_equal(r$stderr, 0)

  expect_warning(r <- brunnerMunzelTest(4:6, 1:3), "variance estimator is zero")
  expect_equal(unname(r$estimate), 0)
  expect_equal(unname(r$statistic), -Inf)

  # requesting the permutation test directly must not warn
  expect_silent(brunnerMunzelTest(1:3, 4:6, method = "permutation"))
})


test_that("completely tied data give a degenerate but finite result", {

  expect_warning(r <- brunnerMunzelTest(rep(2, 5), rep(2, 6)),
                 "variance estimator is zero")
  expect_equal(unname(r$estimate), 0.5)
  expect_equal(unname(r$statistic), 0)
  expect_equal(r$p.value, 1)
  expect_equal(as.vector(r$conf.int), c(NA_real_, NA_real_))
})


test_that("the two-sided permutation p-value uses |T*| >= |T|", {

  # Documented convention, matching the reference implementation of
  # Neubert & Brunner (2007). It is NOT twice the smaller one-sided tail: the
  # permutation distribution here is asymmetric, and 2 * min(tail) would be 0.2.
  expect_warning(r <- brunnerMunzelTest(c(0, 0), c(1, 1, 1)),
                 "variance estimator is zero")
  expect_equal(r$p.value, 0.1)
  expect_equal(brunnerMunzelTest(c(0, 0), c(1, 1, 1),
                                 method = "permutation",
                                 alternative = "greater")$p.value, 0.1)
  expect_equal(brunnerMunzelTest(c(0, 0), c(1, 1, 1),
                                 method = "permutation",
                                 alternative = "less")$p.value, 1)
})


test_that("ordered factors are accepted, unordered factors and text are not", {

  xo <- ordered(c("low", "low", "mid", "high", "high"),
                levels = c("low", "mid", "high"))
  yo <- ordered(c("mid", "high", "high", "high"),
                levels = c("low", "mid", "high"))

  expect_equal(brunnerMunzelTest(xo, yo)$p.value,
               brunnerMunzelTest(c(1, 1, 2, 3, 3), c(2, 3, 3, 3))$p.value)

  expect_error(brunnerMunzelTest(factor(c("a", "b", "a")), 1:3), "ordered factors")
  expect_error(brunnerMunzelTest(c("a", "b", "c"), 1:3), "ordered factors")
  expect_error(brunnerMunzelTest(c(TRUE, FALSE, TRUE), 1:3), "ordered factors")

  # same labels, opposite level order: the integer codes would be mirrored
  yRev <- ordered(as.character(yo), levels = c("high", "mid", "low"))
  expect_error(brunnerMunzelTest(xo, yRev), "identical levels")

  # mixing an ordered factor with a plain numeric vector has no common scale
  expect_error(brunnerMunzelTest(xo, c(1, 2, 3)), "identical levels")
  expect_error(brunnerMunzelTest(c(1, 2, 3), yo), "identical levels")
})


test_that("large separated samples degrade to sampling instead of erroring", {

  # C(80, 40) = 1.1e23 splits, far beyond the enumeration ceiling; a default
  # call must still return rather than demand 'exact = FALSE'
  set.seed(4)
  expect_warning(r <- brunnerMunzelTest(1:40, 41:80), "variance estimator is zero")
  expect_match(r$method, "resamples", fixed = TRUE)
  expect_equal(r$p.value, 1 / (10000 + 1))
  expect_equal(as.vector(r$conf.int), c(NA_real_, NA_real_))

  # small separated samples are still enumerated exactly
  expect_warning(r <- brunnerMunzelTest(1:3, 4:6), "variance estimator is zero")
  expect_match(r$method, "exact studentized permutation", fixed = TRUE)
})


test_that("missing values are removed", {

  expect_equal(brunnerMunzelTest(c(x, NA), c(y, NA, NA))$p.value,
               brunnerMunzelTest(x, y)$p.value)
})


test_that("the confidence interval is clipped to the unit interval", {

  # samples overlap at 5, so the estimator is not degenerate, but the raw Wald
  # upper limit is 1.045 and must be clipped
  a <- c(1, 2, 3, 4, 5)
  b <- c(5, 6, 7, 8, 9)
  ci <- brunnerMunzelTest(a, b)$conf.int
  expect_equal(as.vector(ci), c(0.914776, 1), tolerance = 1e-5)
  expect_gte(min(ci), 0)
  expect_lte(max(ci), 1)

  ci <- brunnerMunzelTest(x, y, conf.level = 0.90)$conf.int
  expect_equal(as.vector(ci), c(0.629098293325597, 0.948823784596481),
               tolerance = 1e-12)
})


test_that("the formula method matches the default method", {

  d <- data.frame(score = c(x, y),
                  grp = rep(c("a", "b"), c(length(x), length(y))))

  f <- brunnerMunzelTest(score ~ grp, data = d)
  g <- brunnerMunzelTest(x, y)

  expect_equal(unname(f$statistic), unname(g$statistic))
  expect_equal(f$p.value, g$p.value)
  expect_equal(f$data.name, "score ~ grp")

  # subset is evaluated in 'data'
  d2 <- rbind(d, data.frame(score = 99, grp = "a"))
  expect_equal(brunnerMunzelTest(score ~ grp, data = d2, subset = score < 99)$p.value,
               g$p.value)

  # more than two groups is not a two-sample design
  d3 <- rbind(d, data.frame(score = 1:3, grp = "c"))
  expect_error(brunnerMunzelTest(score ~ grp, data = d3), "not allowed")
})


test_that("invalid arguments are rejected", {

  expect_error(brunnerMunzelTest(1, y), "at least 2")
  expect_error(brunnerMunzelTest(x, NA_real_), "at least 2")
  expect_error(brunnerMunzelTest(x, y, p0 = 0), "'p0'")
  expect_error(brunnerMunzelTest(x, y, p0 = 1.2), "'p0'")
  expect_error(brunnerMunzelTest(x, y, conf.level = 1), "'conf.level'")
  expect_error(brunnerMunzelTest(x, y, nPerm = 0), "'nPerm'")
  expect_error(brunnerMunzelTest(x, y, nPerm = 1.5), "'nPerm'")
  expect_error(brunnerMunzelTest(x, y, nPerm = Inf), "'nPerm'")
  expect_error(brunnerMunzelTest(x, y, nPerm = NA), "'nPerm'")
  expect_error(brunnerMunzelTest(x, y, nPerm = .Machine$integer.max + 1), "'nPerm'")

  expect_error(brunnerMunzelTest(x, y, exact = NA), "'exact'")
  expect_error(brunnerMunzelTest(x, y, exact = 2), "'exact'")
  expect_error(brunnerMunzelTest(x, y, exact = c(TRUE, FALSE)), "'exact'")
  expect_error(brunnerMunzelTest(x, y, method = "bootstrap"), "should be one of")

  # the permutation null forces p0 = 0.5
  expect_error(brunnerMunzelTest(x, y, p0 = 0.6, method = "permutation"),
               "only for 'p0 = 0.5'", fixed = TRUE)

  # forced enumeration is refused rather than silently started
  set.seed(9)
  big <- rnorm(80)
  expect_error(brunnerMunzelTest(big[1:40], big[41:80], method = "permutation",
                                 exact = TRUE), "exact = FALSE", fixed = TRUE)

  # 'exact' is meaningless without the permutation method
  expect_warning(brunnerMunzelTest(x, y, exact = TRUE), "is ignored")
})
