
# ── Tests for varCI() ────────────────────────────────────────────────────────

library(testthat)

# shared fixtures
x_norm  <- c(14.816, 14.863, 14.814, 14.998, 14.965, 14.824, 14.884, 14.838,
             14.916, 15.021, 14.874, 14.856, 14.860, 14.772, 14.980, 14.919)
x_bonett <- c(15.83, 16.01, 16.24, 16.42, 15.33, 15.44, 16.88, 16.31)
x_nas   <- c(1, 2, NA, 4, 5)

# ── return structure ──────────────────────────────────────────────────────────

test_that("result is a named numeric vector with var/lci/uci", {
  out <- varCI(x_norm)
  expect_type(out, "double")
  expect_named(out, c("var", "lci", "uci"))
})

test_that("var equals sample variance", {
  out <- varCI(x_norm)
  expect_equal(out[["var"]], var(x_norm))
})

test_that("lci <= var <= uci", {
  for (method in c("classic", "bonett")) {
    out <- varCI(x_norm, method = method)
    expect_true(out[["lci"]] <= out[["var"]])
    expect_true(out[["var"]] <= out[["uci"]])
  }
})

# ── classic method ────────────────────────────────────────────────────────────

test_that("classic: CI widens with lower confidence level", {
  ci_95 <- varCI(x_norm, conf.level = 0.95, method = "classic")
  ci_90 <- varCI(x_norm, conf.level = 0.90, method = "classic")
  expect_true(
    (ci_90[["uci"]] - ci_90[["lci"]]) < (ci_95[["uci"]] - ci_95[["lci"]])
  )
})



test_that("classic: CI contains true variance for normal data", {
  set.seed(42)
  x <- rnorm(100, mean = 0, sd = 2)
  out <- varCI(x, conf.level = 0.95, method = "classic")
  expect_true(out[["lci"]] <= 4 && 4 <= out[["uci"]])
})

test_that("classic: known result for x_norm at 90%", {
  out <- varCI(x_norm, conf.level = 0.9, method = "classic")
  # chi-square CI: verify manually
  df <- length(x_norm) - 1
  v  <- var(x_norm)
  expect_equal(out[["lci"]], df * v / qchisq(0.05, df, lower.tail = FALSE),
               tolerance = 1e-10)
  expect_equal(out[["uci"]], df * v / qchisq(0.05, df),
               tolerance = 1e-10)
})

# ── bonett method ─────────────────────────────────────────────────────────────

test_that("bonett: known results from Bonett (2006) paper at 90%", {
  out <- sqrt(varCI(x_bonett, method = "bonett", conf.level = 0.90))
  expect_equal(out[["var"]]^2, var(x_bonett))
  expect_equal(out[["lci"]], 0.3592, tolerance = 1e-3)
  expect_equal(out[["uci"]], 0.9359, tolerance = 1e-3)
})

test_that("bonett: known results from Bonett (2006) paper at 95%", {
  out <- sqrt(varCI(x_bonett, method = "bonett", conf.level = 0.95))
  expect_equal(out[["lci"]], 0.3263, tolerance = 1e-3)
  expect_equal(out[["uci"]], 1.0841, tolerance = 1e-3)
})

test_that("bonett: known results from Bonett (2006) paper at 99%", {
  out <- sqrt(varCI(x_bonett, method = "bonett", conf.level = 0.99))
  expect_equal(out[["lci"]], 0.2607, tolerance = 1e-3)
  expect_equal(out[["uci"]], 1.5109, tolerance = 1e-3)
})

test_that("bonett: stop for n <= 4", {
  expect_error(varCI(1:4, method = "bonett"), regexp = "n > 4")
})

# ── bootstrap method ──────────────────────────────────────────────────────────

test_that("boot: result has correct structure", {
  set.seed(1)
  out <- varCI(x_norm, method = "boot", R = 199)
  expect_named(out, c("var", "lci", "uci"))
  expect_equal(out[["var"]], var(x_norm))
})

test_that("boot: lci <= var <= uci", {
  set.seed(1)
  out <- varCI(x_norm, method = "boot", R = 199)
  expect_true(out[["lci"]] <= out[["var"]])
  expect_true(out[["var"]] <= out[["uci"]])
})

test_that("boot: type = 'norm' produces valid interval", {
  set.seed(1)
  out <- varCI(x_norm, method = "boot", type = "norm", R = 199)
  expect_true(is.finite(out[["lci"]]))
  expect_true(is.finite(out[["uci"]]))
})

test_that("boot: type = 'perc' produces valid interval", {
  set.seed(1)
  out <- varCI(x_norm, method = "boot", type = "perc", R = 199)
  expect_true(is.finite(out[["lci"]]))
  expect_true(is.finite(out[["uci"]]))
})

# ── sides ─────────────────────────────────────────────────────────────────────

test_that("sides = 'left': uci is Inf", {
  out <- varCI(x_norm, sides = "left")
  expect_equal(out[["uci"]], Inf)
  expect_true(is.finite(out[["lci"]]))
})

test_that("sides = 'right': lci is 0", {
  out <- varCI(x_norm, sides = "right")
  expect_equal(out[["lci"]], 0)
  expect_true(is.finite(out[["uci"]]))
})

# Zeile 130-134: fix falsche Richtung
test_that("sides = 'left': lci is higher than two.sided lci", {
  two  <- varCI(x_norm, sides = "two.sided", conf.level = 0.95)
  left <- varCI(x_norm, sides = "left",      conf.level = 0.95)
  expect_true(left[["lci"]] > two[["lci"]])
})


# ── na.rm ─────────────────────────────────────────────────────────────────────

test_that("na.rm = TRUE removes NAs and produces valid result", {
  out <- varCI(x_nas, na.rm = TRUE)
  expect_named(out, c("var", "lci", "uci"))
  expect_equal(out[["var"]], var(c(1, 2, 4, 5)))
})

test_that("na.rm = FALSE with NAs returns NA for var", {
  out <- varCI(x_nas, na.rm = FALSE)
  expect_true(is.na(out[["var"]]))
})

# ── input validation ──────────────────────────────────────────────────────────

test_that("non-numeric input throws error", {
  expect_error(varCI(c("a", "b", "c")), regexp = "numeric")
})

test_that("fewer than 2 observations throws error", {
  expect_error(varCI(42), regexp = "at least two")
})