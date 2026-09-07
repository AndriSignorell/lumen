library(testthat)
library(lumen)

# ===============================================================
# medianCI TESTS
# ===============================================================

# helper: the structural invariants every result has to satisfy.
# Expectations raised inside a helper are counted by testthat, so this
# stays a function rather than being inlined everywhere.
.check_medianCI <- function(res) {

  expect_true(is.numeric(res))
  expect_true(all(c("median", "lci", "uci") %in% names(res)))

  expect_true(is.na(res[["lci"]]) || is.infinite(res[["lci"]]) ||
                res[["lci"]] <= res[["median"]])
  expect_true(is.na(res[["uci"]]) || is.infinite(res[["uci"]]) ||
                res[["uci"]] >= res[["median"]])
  expect_true(is.na(res[["lci"]]) || is.na(res[["uci"]]) ||
                res[["lci"]] <= res[["uci"]])

  invisible(TRUE)
}

set.seed(448)
x.na <- c(rnorm(100), NA)
x    <- x.na[!is.na(x.na)]


test_that("medianCI: the exact method has the documented structure", {

  res <- medianCI(x, method = "exact")

  .check_medianCI(res)
  expect_equal(unname(res[["median"]]), median(x), tolerance = 1e-10)

  # the achieved level is attached as an attribute
  expect_false(is.null(attr(res, "conf.level")))
  expect_gte(attr(res, "conf.level"), 0.90)
  expect_lte(attr(res, "conf.level"), 1.00)
})


test_that("medianCI: the boot method stays close to the exact one", {

  res.ex <- medianCI(x, method = "exact")

  set.seed(1)
  res.boot <- medianCI(x, method = "boot")

  .check_medianCI(res.boot)
  expect_equal(unname(res.boot[["median"]]), median(x), tolerance = 1e-10)
  expect_lt(abs(res.boot[["lci"]] - res.ex[["lci"]]), 0.15)
  expect_lt(abs(res.boot[["uci"]] - res.ex[["uci"]]), 0.15)
})


test_that("medianCI: na.rm = TRUE removes the missing values", {

  res.ex   <- medianCI(x, method = "exact")
  res.narm <- medianCI(x.na, na.rm = TRUE)

  expect_equal(unname(res.narm[["median"]]), median(x), tolerance = 1e-10)
  expect_equal(unname(res.narm[["lci"]]), unname(res.ex[["lci"]]),
               tolerance = 1e-10)
})


test_that("medianCI: a higher conf.level gives a wider interval", {

  res.95 <- medianCI(x, conf.level = 0.95, method = "exact")
  res.99 <- medianCI(x, conf.level = 0.99, method = "exact")

  # the exact interval is built from order statistics, so the achieved
  # level is discrete: both levels can land on the same pair of statistics
  same.level <- attr(res.99, "conf.level") == attr(res.95, "conf.level")

  expect_true(res.99[["lci"]] <= res.95[["lci"]] || same.level)
  expect_true(res.99[["uci"]] >= res.95[["uci"]] || same.level)
})


test_that("medianCI: one-sided intervals open the free side", {

  res.left  <- medianCI(x, sides = "left")
  res.right <- medianCI(x, sides = "right")

  expect_true(is.infinite(res.left[["uci"]]))
  expect_true(is.finite(res.left[["lci"]]))

  expect_true(is.infinite(res.right[["lci"]]))
  expect_true(is.finite(res.right[["uci"]]))
})


test_that("medianCI: fewer than six observations cannot be bounded", {

  res <- medianCI(x = c(1, 2, 3), conf.level = 0.95, method = "exact")

  expect_true(is.infinite(res[["lci"]]) || is.infinite(res[["uci"]]) ||
                attr(res, "conf.level") == 1)
})


test_that("medianCI: symmetric data give a symmetric interval", {

  res <- medianCI(-5:5, method = "exact")

  expect_equal(unname(res[["median"]]), 0, tolerance = 1e-10)
  expect_equal(unname(res[["lci"]]), -unname(res[["uci"]]), tolerance = 1e-10)
})


test_that("medianCI: constant data collapse the interval to a point", {

  res <- medianCI(x = rep(5, 20), method = "exact")

  expect_equal(unname(res[["median"]]), 5, tolerance = 1e-10)
  expect_equal(unname(res[["lci"]]), 5)
  expect_equal(unname(res[["uci"]]), 5)
})


test_that("medianCI: every supported boot type returns a triple", {

  set.seed(42)

  for (btype in c("norm", "basic", "perc", "bca")) {

    res <- medianCI(x, method = "boot", type = btype)

    expect_true(is.numeric(res), info = btype)
    expect_length(res, 3L)
  }
})


test_that("medianCI: an unsupported boot type warns and returns NA limits", {

  expect_warning(res <- medianCI(x, method = "boot", type = "stud"))
  expect_true(anyNA(res[c("lci", "uci")]))
})


test_that("medianCI: the result names never change", {

  for (m in c("exact", "boot")) {
    set.seed(1)
    expect_identical(names(medianCI(x, method = m)),
                     c("median", "lci", "uci"))
  }
})


test_that("medianCI: a single observation is its own median", {

  res <- medianCI(x = 42, method = "exact")
  expect_equal(unname(res[["median"]]), 42)
})
