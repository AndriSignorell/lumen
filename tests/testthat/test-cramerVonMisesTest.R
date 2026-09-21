library(testthat)
library(lumen)

test_that("cramerVonMisesTest: returns htest", {
  set.seed(1)
  res <- cramerVonMisesTest(rnorm(50))
  expect_s3_class(res, "htest")
})

test_that("cramerVonMisesTest: result has statistic and p.value", {
  set.seed(1)
  res <- cramerVonMisesTest(rnorm(50))
  expect_true(!is.null(res$statistic))
  expect_true(!is.null(res$p.value))
})

test_that("cramerVonMisesTest: p.value in [0,1]", {
  set.seed(1)
  res <- cramerVonMisesTest(rnorm(50))
  expect_true(res$p.value >= 0 && res$p.value <= 1)
})

test_that("cramerVonMisesTest: normal data gives large p-value", {
  set.seed(42)
  res <- cramerVonMisesTest(rnorm(200))
  expect_gt(res$p.value, 0.05)
})

test_that("cramerVonMisesTest: uniform data gives small p-value", {
  set.seed(1)
  res <- cramerVonMisesTest(runif(200))
  expect_lt(res$p.value, 0.05)
})

test_that("cramerVonMisesTest: statistic > 0", {
  set.seed(1)
  expect_gt(cramerVonMisesTest(rnorm(30))$statistic, 0)
})

test_that("cramerVonMisesTest: n < 8 throws error", {
  expect_error(cramerVonMisesTest(rnorm(7)))
})

test_that("cramerVonMisesTest: NAs are removed", {
  set.seed(1)
  x <- c(rnorm(50), NA, NA)
  # should not error, NAs silently dropped
  res <- cramerVonMisesTest(x)
  expect_s3_class(res, "htest")
})

test_that("cramerVonMisesTest: method string correct", {
  res <- cramerVonMisesTest(rnorm(20))
  expect_equal(res$method, "Cramer-von Mises normality test")
})


test_that("cramerVonMisesTest: identical to nortest::cvm.test", {
  skip_if_not_installed("nortest")

  set.seed(1)
  for (dd in list(rnorm(100, 5, 3), runif(80), rt(60, 3))) {
    a <- cramerVonMisesTest(dd)
    b <- nortest::cvm.test(dd)
    expect_equal(unname(a$statistic), unname(b$statistic), tolerance = 1e-12)
    expect_equal(a$p.value, b$p.value, tolerance = 1e-12)
  }
})


# -- added --------------------------------------------------------------------

test_that("cramerVonMisesTest: identical to nortest::cvm.test in every branch", {
  skip_if_not_installed("nortest")
  set.seed(21)
  cases <- list(rnorm(15), rnorm(40), rt(50, 4), rexp(30), rexp(80),
                runif(60), rlnorm(40), rchisq(25, 3))
  for (i in 1:30) cases[[length(cases) + 1L]] <- rnorm(10 + i)
  WW <- numeric(0)
  for (i in seq_along(cases)) {
    a <- cramerVonMisesTest(cases[[i]])
    b <- nortest::cvm.test(cases[[i]])
    expect_equal(unname(a$statistic), unname(b$statistic), info = i)
    expect_equal(a$p.value, b$p.value, info = i)
    n <- length(cases[[i]])
    WW <- c(WW, (1 + 0.5 / n) * unname(a$statistic))
  }
  # all four approximation branches below 1.1 were reached
  expect_true(any(WW < 0.0275))
  expect_true(any(WW >= 0.0275 & WW < 0.051))
  expect_true(any(WW >= 0.051 & WW < 0.092))
  expect_true(any(WW >= 0.092 & WW < 1.1))
})

test_that("cramerVonMisesTest: beyond WW = 1.1 the p-value is capped with a warning", {
  set.seed(5)
  x <- c(rep(0, 200), rexp(200, 0.01))
  expect_warning(r <- cramerVonMisesTest(x), "7.37e-10")
  expect_equal(r$p.value, 7.37e-10)
  if (requireNamespace("nortest", quietly = TRUE))
    expect_equal(unname(r$statistic),
                 unname(suppressWarnings(nortest::cvm.test(x))$statistic))
})

test_that("cramerVonMisesTest: invariant under location and scale", {
  set.seed(3)
  x <- rgamma(40, 2)
  expect_equal(cramerVonMisesTest(x)$statistic,
               cramerVonMisesTest(-5 + 0.2 * x)$statistic)
})

test_that("cramerVonMisesTest: input checks", {
  expect_error(cramerVonMisesTest(letters), "numeric")
  expect_error(cramerVonMisesTest(rep(2, 10)), "identical")
  expect_error(cramerVonMisesTest(c(1:10, -Inf)), "infinite")
})
