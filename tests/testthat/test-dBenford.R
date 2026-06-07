library(testthat)
library(lumen)

# --- dbenford ---

test_that("dbenford: first digit probabilities sum to 1", {
  expect_equal(sum(dbenford(1:9)), 1, tolerance = 1e-10)
})

test_that("dbenford: P(D=1) = log10(2) ~ 0.301", {
  expect_equal(dbenford(1), log10(2), tolerance = 1e-10)
})

test_that("dbenford: P(D=9) = log10(10/9)", {
  expect_equal(dbenford(9), log10(10/9), tolerance = 1e-10)
})

test_that("dbenford: out-of-range gives 0", {
  expect_equal(dbenford(0), 0)
  expect_equal(dbenford(10), 0)
})

test_that("dbenford: NA input gives NA", {
  expect_true(is.na(dbenford(NA)))
})

test_that("dbenford: non-integer gives 0", {
  expect_equal(dbenford(1.5), 0)
})

test_that("dbenford log=TRUE returns log of density", {
  expect_equal(dbenford(1, log = TRUE), log(log10(2)), tolerance = 1e-10)
})

test_that("dbenford ndigits=2: probabilities sum to 1", {
  expect_equal(sum(dbenford(10:99, ndigits = 2)), 1, tolerance = 1e-10)
})

test_that("dbenford ndigits=2: out-of-range gives 0", {
  expect_equal(dbenford(9,  ndigits = 2), 0)
  expect_equal(dbenford(100, ndigits = 2), 0)
})

# --- pbenford ---

test_that("pbenford: CDF at 1 equals dbenford(1)", {
  expect_equal(pbenford(1), dbenford(1), tolerance = 1e-10)
})

test_that("pbenford: CDF at 9 equals 1", {
  expect_equal(pbenford(9), 1, tolerance = 1e-10)
})

test_that("pbenford: CDF at 0 equals 0", {
  expect_equal(pbenford(0), 0)
})

test_that("pbenford: non-decreasing", {
  p <- pbenford(1:9)
  expect_true(all(diff(p) >= 0))
})

test_that("pbenford: log.p=TRUE returns log of CDF", {
  expect_equal(pbenford(5, log.p = TRUE), log(pbenford(5)), tolerance = 1e-10)
})

# --- qbenford ---

test_that("qbenford: p=0 returns 1 (lowerlimit)", {
  expect_equal(qbenford(0), 1L)
})

test_that("qbenford: p=1 returns 9 (upperlimit)", {
  expect_equal(qbenford(1), 9L)
})

test_that("qbenford: pbenford(qbenford(p)) == p for exact CDF values", {
  # Benford is discrete: roundtrip exact only at actual CDF values = log10(1+d)
  p_exact <- log10(1 + 1:9)  # = pbenford(1:9) = c(0.301, 0.477, ..., 1)
  expect_equal(pbenford(qbenford(p_exact)), p_exact, tolerance = 1e-10)
})

test_that("qbenford: invalid p throws error", {
  expect_error(qbenford(-0.1))
  expect_error(qbenford(1.1))
})


test_that("dbenford: invalid ndigits throws error", {
  expect_error(dbenford(1, ndigits = 0))
  expect_error(dbenford(1, ndigits = 3))
})

test_that("dbenford: invalid log argument throws error", {
  expect_error(dbenford(1, log = 1))
})

test_that("dbenford: NaN propagates", {
  expect_true(is.nan(dbenford(NaN)))
})

test_that("pbenford ndigits=2 works at bounds", {
  expect_equal(pbenford(9, ndigits = 2), 0)
  expect_equal(pbenford(99, ndigits = 2), 1)
})

test_that("qbenford preserves NA and NaN", {
  expect_true(is.na(qbenford(NA)))
  expect_true(is.na(qbenford(NaN)))
})

test_that("rbenford returns values in support", {
  set.seed(1)
  x <- rbenford(1000)
  expect_true(all(x %in% 1:9))
})

test_that("rbenford ndigits=2 returns values in support", {
  set.seed(1)
  x <- rbenford(1000, ndigits = 2)
  expect_true(all(x %in% 10:99))
})

test_that("rbenford invalid n throws error", {
  expect_error(rbenford(-1))
})

