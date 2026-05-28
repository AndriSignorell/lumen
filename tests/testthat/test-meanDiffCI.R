
# ── Tests for meanDiffCI() ───────────────────────────────────────────────────

library(testthat)

# shared fixtures
x <- mtcars[mtcars$am == 0, "mpg"]   # automatic
y <- mtcars[mtcars$am == 1, "mpg"]   # manual

x_paired <- c(2.1, 3.4, 1.9, 4.2, 3.8)
y_paired <- c(1.8, 3.0, 2.1, 3.9, 3.5)

x_nas <- c(1, 2, NA, 4, 5)
y_nas <- c(2, NA, 4, 5, 6)

# ── return structure ──────────────────────────────────────────────────────────

test_that("result is a named numeric vector with meandiff/lci/uci", {
  out <- meanDiffCI(x, y)
  expect_type(out, "double")
  expect_named(out, c("meandiff", "lci", "uci"))
})

test_that("meandiff equals mean(x) - mean(y)", {
  out <- meanDiffCI(x, y)
  expect_equal(out[["meandiff"]], mean(x) - mean(y))
})

test_that("lci <= meandiff <= uci", {
  out <- meanDiffCI(x, y)
  expect_true(out[["lci"]] <= out[["meandiff"]])
  expect_true(out[["meandiff"]] <= out[["uci"]])
})

# ── classic method ────────────────────────────────────────────────────────────

test_that("classic: matches t.test result", {
  out <- meanDiffCI(x, y, method = "classic")
  tt  <- t.test(x, y)
  expect_equal(out[["lci"]], tt$conf.int[1], tolerance = 1e-10)
  expect_equal(out[["uci"]], tt$conf.int[2], tolerance = 1e-10)
})

test_that("classic: var.equal = TRUE matches pooled t.test", {
  out <- meanDiffCI(x, y, method = "classic", var.equal = TRUE)
  tt  <- t.test(x, y, var.equal = TRUE)
  expect_equal(out[["lci"]], tt$conf.int[1], tolerance = 1e-10)
  expect_equal(out[["uci"]], tt$conf.int[2], tolerance = 1e-10)
})

test_that("classic: CI widens with higher confidence level", {
  ci_95 <- meanDiffCI(x, y, conf.level = 0.95)
  ci_99 <- meanDiffCI(x, y, conf.level = 0.99)
  expect_true(
    (ci_99[["uci"]] - ci_99[["lci"]]) >
      (ci_95[["uci"]] - ci_95[["lci"]])
  )
})

# ── paired ────────────────────────────────────────────────────────────────────

test_that("paired: meandiff equals mean(x - y)", {
  out <- meanDiffCI(x_paired, y_paired, paired = TRUE)
  expect_equal(out[["meandiff"]], mean(x_paired - y_paired))
})

test_that("paired: matches paired t.test result", {
  out <- meanDiffCI(x_paired, y_paired, paired = TRUE)
  tt  <- t.test(x_paired, y_paired, paired = TRUE)
  expect_equal(out[["lci"]], tt$conf.int[1], tolerance = 1e-10)
  expect_equal(out[["uci"]], tt$conf.int[2], tolerance = 1e-10)
})

test_that("paired: unequal lengths throw error", {
  expect_error(
    meanDiffCI(x_paired, c(1, 2, 3), paired = TRUE),
    regexp = "equal length"
  )
})

# ── bootstrap method ──────────────────────────────────────────────────────────

test_that("boot: result has correct structure", {
  set.seed(1)
  out <- meanDiffCI(x, y, method = "boot", R = 199)
  expect_named(out, c("meandiff", "lci", "uci"))
})

test_that("boot: meandiff matches mean(x) - mean(y)", {
  set.seed(1)
  out <- meanDiffCI(x, y, method = "boot", R = 199)
  expect_equal(out[["meandiff"]], mean(x) - mean(y), tolerance = 1e-10)
})

test_that("boot: lci <= meandiff <= uci", {
  set.seed(1)
  out <- meanDiffCI(x, y, method = "boot", R = 199)
  expect_true(out[["lci"]] <= out[["meandiff"]])
  expect_true(out[["meandiff"]] <= out[["uci"]])
})

test_that("boot: type = 'norm' produces finite interval", {
  set.seed(1)
  out <- meanDiffCI(x, y, method = "boot", type = "norm", R = 199)
  expect_true(is.finite(out[["lci"]]))
  expect_true(is.finite(out[["uci"]]))
})

test_that("boot: type = 'perc' produces finite interval", {
  set.seed(1)
  out <- meanDiffCI(x, y, method = "boot", type = "perc", R = 199)
  expect_true(is.finite(out[["lci"]]))
  expect_true(is.finite(out[["uci"]]))
})

test_that("boot paired: meandiff equals mean(x - y)", {
  set.seed(1)
  out <- meanDiffCI(x_paired, y_paired, method = "boot", paired = TRUE, R = 199)
  expect_equal(out[["meandiff"]], mean(x_paired - y_paired), tolerance = 1e-10)
})

# ── sides ─────────────────────────────────────────────────────────────────────

test_that("sides = 'left': uci is Inf", {
  out <- meanDiffCI(x, y, sides = "left")
  expect_equal(out[["uci"]], Inf)
  expect_true(is.finite(out[["lci"]]))
})

test_that("sides = 'right': lci is -Inf", {
  out <- meanDiffCI(x, y, sides = "right")
  expect_equal(out[["lci"]], -Inf)
  expect_true(is.finite(out[["uci"]]))
})

test_that("sides = 'left': lci is higher than two.sided lci", {
  two  <- meanDiffCI(x, y, sides = "two.sided", conf.level = 0.95)
  left <- meanDiffCI(x, y, sides = "left",      conf.level = 0.95)
  expect_true(left[["lci"]] > two[["lci"]])
})

# ── na.rm ─────────────────────────────────────────────────────────────────────

test_that("na.rm = TRUE removes NAs and produces valid result", {
  out <- meanDiffCI(x_nas, y_nas, na.rm = TRUE)
  expect_named(out, c("meandiff", "lci", "uci"))
  expect_equal(out[["meandiff"]], mean(c(1, 2, 4, 5)) - mean(c(2, 4, 5, 6)))
})

test_that("na.rm = FALSE with NAs returns NA for meandiff", {
  out <- meanDiffCI(x_nas, y_nas, na.rm = FALSE)
  expect_true(is.na(out[["meandiff"]]))
})

# ── input validation ──────────────────────────────────────────────────────────

test_that("non-numeric x throws error", {
  expect_error(meanDiffCI(c("a", "b", "c"), y), regexp = "numeric")
})

test_that("non-numeric y throws error", {
  expect_error(meanDiffCI(x, c("a", "b", "c")), regexp = "numeric")
})

test_that("fewer than 2 observations in x throws error", {
  expect_error(meanDiffCI(1, y), regexp = "at least two")
})

test_that("fewer than 2 observations in y throws error", {
  expect_error(meanDiffCI(x, 1), regexp = "at least two")
})

