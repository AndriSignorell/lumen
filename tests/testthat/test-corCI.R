test_that("corCI() follows the package argument conventions", {

  a <- names(formals(corCI))

  expect_false("alternative" %in% a)
  expect_true("sides" %in% a)

  # data, then the inference block (design_rules 4.1)
  expect_equal(a, c("rho", "n", "conf.level", "sides"))

  # a dedicated CI function keeps 0.95 as its default
  expect_equal(formals(corCI)$conf.level, 0.95)
})


test_that("the numbers are unchanged by the rename", {

  # reference values computed from Fisher's z by hand
  ci <- corCI(0.5, n = 30)
  expect_named(ci, c("est", "lci", "uci"))
  expect_equal(unname(ci[["est"]]), 0.5)
  expect_equal(unname(ci[["lci"]]), 0.1704313651, tolerance = 1e-8)
  expect_equal(unname(ci[["uci"]]), 0.7289585564, tolerance = 1e-8)

  left <- corCI(0.5, n = 30, sides = "left")
  expect_equal(unname(left[["lci"]]), 0.2286399420, tolerance = 1e-8)

  right <- corCI(0.5, n = 30, sides = "right")
  expect_equal(unname(right[["uci"]]), 0.6992637581, tolerance = 1e-8)
})


test_that("sides names the side carrying the finite bound", {

  two   <- corCI(0.5, n = 30)
  left  <- corCI(0.5, n = 30, sides = "left")
  right <- corCI(0.5, n = 30, sides = "right")

  # a correlation lives in [-1, 1], so the open side stops there
  expect_equal(left[["uci"]], 1)
  expect_equal(right[["lci"]], -1)

  expect_equal(left[["est"]], two[["est"]])

  # a one-sided bound carries the whole alpha, so it is tighter
  expect_gt(left[["lci"]],  two[["lci"]])
  expect_lt(right[["uci"]], two[["uci"]])

  # "left" is the former alternative = "greater"
  expect_equal(left[["lci"]], corCI(0.5, n = 30, conf.level = 0.90)[["lci"]])
  expect_equal(right[["uci"]], corCI(0.5, n = 30, conf.level = 0.90)[["uci"]])
})


test_that("the interval brackets rho and stays inside the range", {

  for (r in seq(-0.9, 0.9, by = 0.3)) {
    ci <- corCI(r, n = 25)
    expect_true(ci[["lci"]] <= r && r <= ci[["uci"]])
    expect_true(ci[["lci"]] >= -1 && ci[["uci"]] <= 1)
  }
})


test_that("a perfect correlation gets NA bounds, not a degenerate interval", {

  # (rho, rho) would rule out every value below 1, which no sample supports
  expect_warning(ci <- corCI(1, n = 30), "perfect correlation")
  expect_equal(unname(ci[["est"]]), 1)
  expect_true(is.na(ci[["lci"]]))
  expect_true(is.na(ci[["uci"]]))

  expect_warning(cim <- corCI(-1, n = 30), "perfect correlation")
  expect_equal(unname(cim[["est"]]), -1)
  expect_true(is.na(cim[["lci"]]))
})


test_that("n = 3 gets NA bounds rather than the whole range", {

  # 1/sqrt(n-3) is infinite there; the old code returned (-1, 1), which
  # looks like a result and is merely the parameter space
  expect_warning(ci <- corCI(0.5, n = 3), "3 observations")
  expect_true(is.na(ci[["lci"]]))
  expect_equal(unname(ci[["est"]]), 0.5)

  expect_silent(corCI(0.5, n = 4))
})


test_that("the arguments are validated", {

  expect_error(corCI(0.5, n = 2), "at least 3")
  expect_error(corCI(0.5, n = "a"), "'n'")
  expect_error(corCI(0.5, n = 30.5), "whole number")
  expect_error(corCI(0.5, n = c(10, 20)), "'n'")

  expect_error(corCI(1.5, n = 30), "between -1 and 1")
  expect_error(corCI(c(0.1, 0.2), n = 30), "single finite")
  expect_error(corCI(NA, n = 30), "single finite")

  expect_error(corCI(0.5, n = 30, conf.level = c(0.9, 0.95)), "conf.level")
  expect_error(corCI(0.5, n = 30, conf.level = 0), "conf.level")
  expect_error(corCI(0.5, n = 30, conf.level = NA), "conf.level")

  expect_error(corCI(0.5, n = 30, conf.level = 0.4, sides = "left"), "0.5")
  expect_silent(corCI(0.5, n = 30, conf.level = 0.4))

  expect_error(corCI(0.5, n = 30, sides = "links"), "two.sided")
})


test_that("corCI() agrees with cor.test() on the same data", {

  set.seed(1)
  x <- rnorm(40)
  y <- 0.6 * x + rnorm(40, sd = 0.8)

  ct <- cor.test(x, y)
  ci <- corCI(unname(ct$estimate), n = 40)

  expect_equal(unname(ci[["lci"]]), ct$conf.int[1], tolerance = 1e-8)
  expect_equal(unname(ci[["uci"]]), ct$conf.int[2], tolerance = 1e-8)

  # and the one-sided correspondence: sides = "left" is greater
  ctg <- cor.test(x, y, alternative = "greater")
  expect_equal(unname(corCI(unname(ct$estimate), n = 40,
                            sides = "left")[["lci"]]),
               ctg$conf.int[1], tolerance = 1e-8)
})
