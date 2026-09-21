library(testthat)
library(lumen)

# Fixed test data: iid errors -> no serial correlation
set.seed(42)
x  <- rep(c(1, -1), 50)
y1 <- 1 + x + rnorm(100)
# AR(1) errors: strong autocorrelation
y2 <- stats::filter(y1, 0.9, method = "recursive")

test_that("breuschGodfreyTest: returns htest", {
  res <- breuschGodfreyTest(y1 ~ x)
  expect_s3_class(res, "htest")
})

test_that("breuschGodfreyTest: result components present", {
  res <- breuschGodfreyTest(y1 ~ x)
  expect_false(is.null(res$statistic))
  expect_false(is.null(res$p.value))
  expect_false(is.null(res$parameter))
})

test_that("breuschGodfreyTest: p.value in [0,1]", {
  res <- breuschGodfreyTest(y1 ~ x)
  expect_true(res$p.value >= 0 && res$p.value <= 1)
})

test_that("breuschGodfreyTest: iid residuals give large p-value", {
  res <- breuschGodfreyTest(y1 ~ x, order = 1)
  expect_gt(res$p.value, 0.05)
})

test_that("breuschGodfreyTest: AR(1) residuals give small p-value", {
  res <- breuschGodfreyTest(y2 ~ x, order = 1)
  expect_lt(res$p.value, 0.05)
})

test_that("breuschGodfreyTest: order=4 gives df=4", {
  res <- breuschGodfreyTest(y1 ~ x, order = 4)
  expect_equal(unname(res$parameter), 4L)
})

test_that("breuschGodfreyTest: type='F' returns F statistic", {
  res <- breuschGodfreyTest(y1 ~ x, type = "f")
  expect_named(res$parameter, c("df1", "df2"))
})

test_that("breuschGodfreyTest: accepts lm object", {
  fit <- lm(y1 ~ x)
  res <- breuschGodfreyTest(fit)
  expect_s3_class(res, "htest")
})

test_that("breuschGodfreyTest: coefficients and vcov present", {
  res <- breuschGodfreyTest(y1 ~ x)
  expect_false(is.null(res$coefficients))
  expect_false(is.null(res$vcov))
})

test_that("breuschGodfreyTest: statistic >= 0", {
  res <- breuschGodfreyTest(y1 ~ x)
  expect_gte(unname(res$statistic), 0)
})


test_that("breuschGodfreyTest: identical to lmtest::bgtest", {
  skip_if_not_installed("lmtest")

  for (ord in c(1L, 4L)) {
    a <- breuschGodfreyTest(y1 ~ x, order = ord)
    b <- lmtest::bgtest(y1 ~ x, order = ord)
    expect_equal(unname(a$statistic), unname(b$statistic), tolerance = 1e-12)
    expect_equal(a$p.value, b$p.value, tolerance = 1e-12)
    expect_equal(a$coefficients, b$coefficients, tolerance = 1e-12)
    expect_equal(a$vcov, b$vcov, tolerance = 1e-12)
  }

  a <- breuschGodfreyTest(y1 ~ x, order = 2, type = "f")
  b <- lmtest::bgtest(y1 ~ x, order = 2, type = "F")
  expect_equal(unname(a$statistic), unname(b$statistic), tolerance = 1e-12)
  expect_equal(a$p.value, b$p.value, tolerance = 1e-12)
})

test_that("breuschGodfreyTest: type is case-insensitive", {
  a <- breuschGodfreyTest(y1 ~ x, type = "F")
  b <- breuschGodfreyTest(y1 ~ x, type = "f")
  expect_equal(a$statistic, b$statistic)
})

test_that("breuschGodfreyTest: invalid order throws error", {
  expect_error(breuschGodfreyTest(y1 ~ x, order = 0), "positive integer")
  expect_error(breuschGodfreyTest(y1 ~ x, order = -1), "positive integer")
})

test_that("breuschGodfreyTest: vcov and df.residual methods work (coeftest)", {
  skip_if_not_installed("lmtest")

  res <- breuschGodfreyTest(y1 ~ x, order = 2)
  expect_equal(vcov(res), res$vcov)

  # the residual df of the auxiliary regression do not depend on which
  # statistic was requested: n - k - order = 100 - 2 - 2
  expect_identical(df.residual(res), res$df.residual)
  expect_equal(df.residual(res), 96)

  res_f <- breuschGodfreyTest(y1 ~ x, order = 2, type = "f")
  expect_equal(df.residual(res_f), unname(res_f$parameter["df2"]))
  expect_identical(df.residual(res), df.residual(res_f))

  ct <- lmtest::coeftest(res)
  expect_equal(nrow(ct), length(res$coefficients))

  # deliberate deviation from lmtest: bgtest() hands back NULL for the
  # chi-squared version, which sends coeftest() to the normal distribution
  expect_equal(colnames(ct)[4L], "Pr(>|t|)")
  expect_equal(colnames(lmtest::coeftest(lmtest::bgtest(y1 ~ x, order = 2)))[4L],
               "Pr(>|z|)")
})


# -- added --------------------------------------------------------------------

set.seed(3)
bb <- data.frame(x = rnorm(50), z = runif(50))
bb$y <- 1 + bb$x + as.numeric(stats::filter(rnorm(50), 0.4, "recursive"))

test_that("identical to lmtest::bgtest with fill = NA and with orderBy", {
  skip_if_not_installed("lmtest")
  for (tp in c("Chisq", "F")) {
    a <- breuschGodfreyTest(y ~ x, data = bb, order = 3, fill = NA, type = tp)
    b <- lmtest::bgtest(y ~ x, data = bb, order = 3, fill = NA, type = tp)
    expect_equal(unname(a$statistic), unname(b$statistic), info = tp)
    expect_equal(a$p.value, b$p.value, info = tp)
    expect_equal(unname(a$parameter), unname(b$parameter), info = tp)

    a <- breuschGodfreyTest(y ~ x, data = bb, order = 2, orderBy = ~ z, type = tp)
    b <- lmtest::bgtest(y ~ x, data = bb, order = 2, order.by = ~ z, type = tp)
    expect_equal(unname(a$statistic), unname(b$statistic), info = tp)
  }
})

test_that("fill = NA drops the first 'order' observations", {
  a <- breuschGodfreyTest(y ~ x, data = bb, order = 3, fill = NA, type = "F")
  # k = 2 regressors + 3 lags, n = 50 - 3
  expect_equal(a$df.residual, 50 - 3 - 2 - 3)
})

test_that("lm input equals the formula, with and without stored x/y", {
  f <- breuschGodfreyTest(y ~ x, data = bb, order = 2)
  l1 <- breuschGodfreyTest(lm(y ~ x, data = bb), order = 2)
  l2 <- breuschGodfreyTest(lm(y ~ x, data = bb, x = TRUE, y = TRUE), order = 2)
  expect_equal(l1$statistic, f$statistic)
  expect_equal(l2$statistic, f$statistic)
  expect_match(l1$data.name, "lm")
})

test_that("lm input: orderBy aligned after subset", {
  fit <- lm(y ~ x, data = bb, subset = z > 0.3)
  a <- breuschGodfreyTest(fit, data = bb, orderBy = ~ z, order = 2)
  b <- breuschGodfreyTest(y ~ x, data = bb[bb$z > 0.3, ], orderBy = ~ z, order = 2)
  expect_equal(a$statistic, b$statistic)
})

test_that("formula subset", {
  a <- breuschGodfreyTest(y ~ x, data = bb, subset = z > 0.3, order = 2)
  b <- breuschGodfreyTest(y ~ x, data = bb[bb$z > 0.3, ], order = 2)
  expect_equal(a$statistic, b$statistic)
})

test_that("coefficients and vcov are labelled", {
  r <- breuschGodfreyTest(y ~ x, data = bb, order = 2)
  nm <- c("(Intercept)", "x", "lag(resid)_1", "lag(resid)_2")
  expect_named(r$coefficients, nm)
  expect_identical(dimnames(vcov(r)), list(nm, nm))
  expect_identical(df.residual(r), r$df.residual)
  expect_identical(r$method,
                   "Breusch-Godfrey test for serial correlation of order up to 2")
})

test_that("argument checks", {
  expect_error(breuschGodfreyTest(y ~ x, data = bb, order = NA), "positive integer")
  expect_error(breuschGodfreyTest(y ~ x, data = bb, order = c(1, 2)), "positive integer")
  expect_error(breuschGodfreyTest(y ~ x, data = bb, order = 50), "smaller than")
  expect_error(breuschGodfreyTest(y ~ x, data = bb, fill = "a"), "'fill'")
  expect_error(breuschGodfreyTest(y ~ x, data = bb, fill = c(0, 0)), "'fill'")
  expect_error(breuschGodfreyTest(y ~ x, data = bb, type = "t"))
  expect_error(breuschGodfreyTest(cbind(y, z) ~ x, data = bb), "single vector")
})

test_that("rank deficiency and too few observations", {
  d2 <- bb; d2$x2 <- 2 * d2$x
  expect_error(breuschGodfreyTest(y ~ x + x2, data = d2), "model matrix is rank deficient")
  expect_error(breuschGodfreyTest(y ~ x, data = bb[1:5, ], order = 3),
               "not enough observations")
})
