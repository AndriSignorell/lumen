library(testthat)
library(lumen)

# Mehta et al. (2003) example
tab_mehta <- matrix(c(7, 12, 8, 3), nrow = 2,
                    dimnames = list(treat = c("vaccine","placebo"),
                                    infection = c("yes","no")))

# Small balanced table
tab_small <- matrix(c(8, 14, 1, 3), nrow = 2)

test_that("barnardTest: returns htest", {
  res <- barnardTest(tab_small)
  expect_s3_class(res, "htest")
})

test_that("barnardTest: p.value in [0,1]", {
  res <- barnardTest(tab_small)
  expect_true(res$p.value >= 0 && res$p.value <= 1)
})

test_that("barnardTest: two-sided p >= one-sided p", {
  res_two  <- barnardTest(tab_mehta, alternative = "two.sided")
  res_less <- barnardTest(tab_mehta, alternative = "less")
  expect_gte(res_two$p.value, res_less$p.value - 1e-6)
})

test_that("barnardTest: non-2x2 matrix throws error", {
  tab3x2 <- matrix(1:6, nrow = 3)
  expect_error(barnardTest(tab3x2))
})

test_that("barnardTest: Mehta example one-sided p < 0.05", {
  res <- barnardTest(tab_mehta, alternative = "less")
  expect_lt(res$p.value, 0.05)
})

test_that("barnardTest: result has method field", {
  res <- barnardTest(tab_small)
  expect_false(is.null(res$method))
})

test_that("barnardTest: fixed=2 (column margins) works", {
  res <- barnardTest(tab_small, fixed = 2)
  expect_s3_class(res, "htest")
})


# -- added --------------------------------------------------------------------

skip_if_not_installed("Exact")

ex <- function(tab, ...)
  Exact::exact.test(tab, to.plot = FALSE, useStoredCSM = FALSE, ...)

test_that("barnardTest: every method equals Exact::exact.test (rows fixed)", {
  for (m in c("csm", "z-pooled", "z-unpooled", "boschloo")) {
    r <- barnardTest(tab_mehta, method = m)
    e <- ex(tab_mehta, method = m, model = "binomial", cond.row = TRUE)
    expect_equal(r$p.value, e$p.value, info = m)
    expect_equal(r$statistic, e$statistic, info = m)
  }
  r <- barnardTest(tab_mehta, method = "santner-snell")
  e <- ex(tab_mehta, method = "santner and snell", model = "binomial",
          cond.row = TRUE)
  expect_equal(r$p.value, e$p.value)
})

test_that("barnardTest: fixed = 2 conditions on the columns", {
  r <- barnardTest(tab_mehta, fixed = 2, method = "z-pooled")
  e <- ex(tab_mehta, method = "z-pooled", model = "binomial", cond.row = FALSE)
  expect_equal(r$p.value, e$p.value)
  # the same as fixing the rows of the transposed table
  expect_equal(r$p.value,
               barnardTest(t(tab_mehta), fixed = 1, method = "z-pooled")$p.value)
})

test_that("barnardTest: fixed = NA uses the multinomial model", {
  r <- barnardTest(tab_small, fixed = NA, method = "z-pooled")
  e <- ex(tab_small, method = "z-pooled", model = "multinomial")
  expect_equal(r$p.value, e$p.value)
  expect_s3_class(r, "htest")
})

test_that("barnardTest: user values in '...' override the defaults", {
  r <- barnardTest(tab_mehta, method = "z-pooled", model = "multinomial")
  e <- ex(tab_mehta, method = "z-pooled", model = "multinomial")
  expect_equal(r$p.value, e$p.value)
  r <- barnardTest(tab_mehta, method = "z-pooled", cond.row = FALSE)
  expect_equal(r$p.value, barnardTest(tab_mehta, method = "z-pooled",
                                      fixed = 2)$p.value)
})

test_that("barnardTest: alternatives", {
  l <- barnardTest(tab_mehta, alternative = "less", method = "z-pooled")
  g <- barnardTest(tab_mehta, alternative = "greater", method = "z-pooled")
  e <- ex(tab_mehta, method = "z-pooled", model = "binomial", cond.row = TRUE,
          alternative = "greater")
  expect_equal(g$p.value, e$p.value)
  expect_identical(l$alternative, "less")
  expect_error(barnardTest(tab_mehta, alternative = "foo"))
})

test_that("barnardTest: two factors give the same result as their table", {
  df <- as.data.frame(as.table(tab_mehta))
  tr <- rep(df$treat, df$Freq)
  inf <- rep(df$infection, df$Freq)
  r <- barnardTest(tr, inf, method = "z-pooled")
  expect_equal(r$p.value, barnardTest(tab_mehta, method = "z-pooled")$p.value)
  expect_identical(r$data.name, "tr and inf")
})

test_that("barnardTest: data.name is a single string", {
  r <- barnardTest(matrix(c(7, 12, 8, 3), nrow = 2, dimnames = list(
    treatment_group = c("vaccine", "placebo"),
    infection_status = c("yes", "no"))), method = "z-pooled")
  expect_length(r$data.name, 1L)
})

test_that("barnardTest: 'fixed' is validated", {
  expect_error(barnardTest(tab_small, fixed = c(1, 2)), "fisher.test")
  expect_error(barnardTest(tab_small, fixed = c(2, 1)), "fisher.test")
  for (bad in list(3, 0, "1", c(1, 1)))
    expect_error(barnardTest(tab_small, fixed = bad), "'fixed' must be",
                 info = format(bad))
  expect_error(barnardTest(tab_small, method = "fisher"))
})
