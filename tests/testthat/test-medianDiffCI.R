# Tests for medianDiffCI()
#
# Reference values for the Bonett-Price interval were computed
# independently in Python (numpy/scipy) from the same order-statistic
# formula: mtcars$mpg by am, n = 19 / 13, c = 5 / 3.

x <- mtcars$mpg[mtcars$am == 0]
y <- mtcars$mpg[mtcars$am == 1]


# -- structure ----------------------------------------------------------------

test_that("medianDiffCI: returns named numeric vector of length 3", {
  r <- medianDiffCI(x, y)
  expect_type(r, "double")
  expect_named(r, c("mediandiff", "lci", "uci"))
  expect_equal(unname(r["mediandiff"]), median(x) - median(y))
})


# -- Bonett-Price reference values -------------------------------------------

test_that("medianDiffCI: bonett-price matches reference values", {
  expect_equal(
    unname(medianDiffCI(x, y)),
    c(-5.5, -10.465132862361095, -0.5348671376389049), tolerance = 1e-10)
  expect_equal(
    unname(medianDiffCI(x, y, conf.level = 0.90)),
    c(-5.5, -9.666870851388186, -1.333129148611814), tolerance = 1e-10)
  expect_equal(
    unname(medianDiffCI(x, y, conf.level = 0.99)),
    c(-5.5, -12.02529068072882, 1.0252906807288191), tolerance = 1e-10)
})


# local cross-check, statpsych not in Suggests

# test_that("medianDiffCI: bonett-price agrees with statpsych::ci.median2", {
#   skip_if_not_installed("statpsych")
#   ref <- statpsych::ci.median2(0.05, x, y)
#   r   <- medianDiffCI(x, y)
#   expect_equal(unname(r["lci"]), unname(ref[1, "LL"]), tolerance = 1e-8)
#   expect_equal(unname(r["uci"]), unname(ref[1, "UL"]), tolerance = 1e-8)
# })


test_that("medianDiffCI: c is clamped to 1 for tiny samples", {
  # n = 2: round(1 - sqrt(2)) = 0 -> must be clamped, result finite
  r <- medianDiffCI(c(1, 3), c(2, 5))
  expect_true(all(is.finite(r)))
  expect_lt(r["lci"], r["uci"])
})


# -- invariances --------------------------------------------------------------

test_that("medianDiffCI: swapping x and y mirrors the interval", {
  a <- medianDiffCI(x, y)
  b <- medianDiffCI(y, x)
  expect_equal(unname(b), unname(c(-a["mediandiff"], -a["uci"], -a["lci"])))
})

test_that("medianDiffCI: shift and scale equivariance", {
  a <- medianDiffCI(x, y)
  expect_equal(unname(medianDiffCI(x + 10, y)), unname(a + 10))
  expect_equal(unname(medianDiffCI(3 * x, 3 * y)), unname(3 * a))
})

test_that("medianDiffCI: higher conf.level gives a wider, nested interval", {
  r90 <- medianDiffCI(x, y, conf.level = 0.90)
  r95 <- medianDiffCI(x, y, conf.level = 0.95)
  r99 <- medianDiffCI(x, y, conf.level = 0.99)
  expect_true(r99["lci"] < r95["lci"] && r95["lci"] < r90["lci"])
  expect_true(r90["uci"] < r95["uci"] && r95["uci"] < r99["uci"])
})


# -- sides --------------------------------------------------------------------

test_that("medianDiffCI: one-sided bounds equal two-sided at 1 - 2 alpha", {
  two <- medianDiffCI(x, y, conf.level = 0.90)

  l <- medianDiffCI(x, y, sides = "left")
  expect_equal(unname(l["lci"]), unname(two["lci"]))
  expect_identical(unname(l["uci"]), Inf)

  r <- medianDiffCI(x, y, sides = "right")
  expect_equal(unname(r["uci"]), unname(two["uci"]))
  expect_identical(unname(r["lci"]), -Inf)
})

test_that("medianDiffCI: sides accepts initial letter", {
  expect_equal(medianDiffCI(x, y, sides = "l"),
               medianDiffCI(x, y, sides = "left"))
})


# -- missing values -----------------------------------------------------------

test_that("medianDiffCI: NA without na.rm returns NA, no error", {
  r <- medianDiffCI(c(x, NA), y)
  expect_named(r, c("mediandiff", "lci", "uci"))
  expect_true(all(is.na(r)))
})

test_that("medianDiffCI: na.rm = TRUE removes NA per group", {
  expect_equal(medianDiffCI(c(x, NA), c(NA, y), na.rm = TRUE),
               medianDiffCI(x, y))
})

test_that("medianDiffCI: na.rm = TRUE removes incomplete pairs", {
  a <- c(1, 4, NA, 6, 8, 3)
  b <- c(2, 2, 5, NA, 7, 1)
  ok <- complete.cases(a, b)
  set.seed(11)
  r1 <- medianDiffCI(a, b, method = "boot", paired = TRUE, na.rm = TRUE,
                     R = 199, type = "perc")
  set.seed(11)
  r2 <- medianDiffCI(a[ok], b[ok], method = "boot", paired = TRUE,
                     R = 199, type = "perc")
  expect_equal(r1, r2)
})


# -- argument validation ------------------------------------------------------

test_that("medianDiffCI: invalid input raises errors", {
  expect_error(medianDiffCI(letters, y), "must be numeric")
  expect_error(medianDiffCI(1, y), "at least two")
  expect_error(medianDiffCI(x, y[1:5], paired = TRUE), "equal length")
  expect_error(medianDiffCI(x[1:13], y, paired = TRUE),
               "method = \"boot\"")
  expect_error(medianDiffCI(x, y, method = "foo"))
  expect_error(medianDiffCI(x, y, conf.level = 1.5))
  expect_error(medianDiffCI(x, y, method = "boot", type = "stud"),
               "not supported")
})


# -- bootstrap ----------------------------------------------------------------

test_that("medianDiffCI: bootstrap is reproducible and brackets estimate", {
  set.seed(1)
  a <- medianDiffCI(x, y, method = "boot", R = 499, type = "perc")
  set.seed(1)
  b <- medianDiffCI(x, y, method = "boot", R = 499, type = "perc")
  expect_identical(a, b)
  expect_equal(unname(a["mediandiff"]), median(x) - median(y))
  expect_true(a["lci"] <= a["mediandiff"] && a["mediandiff"] <= a["uci"])
})

test_that("medianDiffCI: all supported bootstrap types return finite bounds", {
  for (tp in c("perc", "basic", "norm", "bca")) {
    set.seed(3)
    # bca on medians may warn about extreme order statistics
    r <- suppressWarnings(
      medianDiffCI(x, y, method = "boot", R = 999, type = tp))
    expect_true(all(is.finite(r)), info = tp)
    expect_lt(r["lci"], r["uci"])
  }
})

test_that("medianDiffCI: bootstrap resamples within groups", {
  # x is constant: with stratified resampling median(x*) is always 100,
  # so every replicate lies in [100 - max(y), 100 - min(y)]. Mixing the
  # groups would pull median(x*) towards 0 in many replicates.
  set.seed(4)
  yy <- rnorm(15)
  r  <- medianDiffCI(rep(100, 15), yy, method = "boot", R = 499,
                     type = "perc")
  expect_gte(r["lci"], 100 - max(yy))
  expect_lte(r["uci"], 100 - min(yy))
})

test_that("medianDiffCI: paired estimate is difference of medians", {
  d2 <- sleep$extra[sleep$group == 2]
  d1 <- sleep$extra[sleep$group == 1]
  set.seed(5)
  r <- medianDiffCI(d2, d1, method = "boot", paired = TRUE, R = 499,
                    type = "perc")
  expect_equal(unname(r["mediandiff"]), median(d2) - median(d1))
  expect_false(isTRUE(all.equal(unname(r["mediandiff"]),
                                median(d2 - d1))))
})

test_that("medianDiffCI: sides works with bootstrap", {
  set.seed(6)
  r <- medianDiffCI(x, y, method = "boot", sides = "left", R = 499,
                    type = "perc")
  expect_identical(unname(r["uci"]), Inf)
  expect_true(is.finite(r["lci"]))
})


# -- coverage (slow) ----------------------------------------------------------

test_that("medianDiffCI: bonett-price coverage close to nominal", {
  skip_on_cran()
  set.seed(2026)
  # different shapes and scales, true median difference log(2) - 2 log(2)
  hit <- replicate(2000, {
    r <- medianDiffCI(rexp(30, 1), rexp(50, 0.5))
    r["lci"] <= -log(2) && -log(2) <= r["uci"]
  })
  expect_gt(mean(hit), 0.93)
  expect_lt(mean(hit), 0.97)
})
