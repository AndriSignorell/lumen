library(testthat)
library(lumen)

set.seed(42)
x <- rnorm(200, mean = 5, sd = 2)
y <- rnorm(150, mean = 5, sd = 3)

# -----------------------------------------------------------------------
# median_boot_cpp
# -----------------------------------------------------------------------







# -----------------------------------------------------------------------
# mad_boot_cpp
# -----------------------------------------------------------------------

test_that("mad_boot_cpp: returns named numeric vector of length 3", {
  res <- mad_boot_cpp(x, R = 500, alpha = 0.05, seed = 1)
  expect_true(is.numeric(res))
  expect_length(res, 3L)
  expect_named(res, c("est", "lci", "uci"))
})

test_that("mad_boot_cpp: lci <= est <= uci", {
  res <- mad_boot_cpp(x, R = 500, alpha = 0.05, seed = 1)
  expect_lte(res["lci"], res["est"])
  expect_lte(res["est"], res["uci"])
})

test_that("mad_boot_cpp: est matches mad(x)", {
  res <- mad_boot_cpp(x, R = 500, alpha = 0.05, seed = 1)
  expect_equal(unname(res["est"]), mad(x), tolerance = 1e-10)
})

test_that("mad_boot_cpp: custom constant is respected", {
  res1 <- mad_boot_cpp(x, R = 500, constant = 1.4826, seed = 1)
  res2 <- mad_boot_cpp(x, R = 500, constant = 1.0,    seed = 1)
  expect_false(isTRUE(all.equal(res1["est"], res2["est"])))
})

test_that("mad_boot_cpp: seed ensures reproducibility", {
  r1 <- mad_boot_cpp(x, R = 500, seed = 7)
  r2 <- mad_boot_cpp(x, R = 500, seed = 7)
  expect_equal(r1, r2)
})

# -----------------------------------------------------------------------
# mad_diff_boot_cpp
# -----------------------------------------------------------------------

test_that("mad_diff_boot_cpp: returns named numeric vector of length 3", {
  res <- mad_diff_boot_cpp(x, y, R = 500, alpha = 0.05, seed = 1)
  expect_true(is.numeric(res))
  expect_length(res, 3L)
  expect_named(res, c("est", "lci", "uci"))
})

test_that("mad_diff_boot_cpp: est matches mad(x) - mad(y)", {
  res <- mad_diff_boot_cpp(x, y, R = 500, seed = 1)
  expect_equal(unname(res["est"]), mad(x) - mad(y), tolerance = 1e-10)
})

test_that("mad_diff_boot_cpp: lci <= est <= uci", {
  res <- mad_diff_boot_cpp(x, y, R = 500, seed = 1)
  expect_lte(res["lci"], res["est"])
  expect_lte(res["est"], res["uci"])
})

test_that("mad_diff_boot_cpp: identical samples give est near 0", {
  res <- mad_diff_boot_cpp(x, x, R = 500, seed = 1)
  expect_equal(unname(res["est"]), 0, tolerance = 1e-10)
})

# -----------------------------------------------------------------------
# mad_ratio_boot_cpp
# -----------------------------------------------------------------------

test_that("mad_ratio_boot_cpp: returns named numeric vector of length 3", {
  res <- mad_ratio_boot_cpp(x, y, R = 500, alpha = 0.05, seed = 1)
  expect_true(is.numeric(res))
  expect_length(res, 3L)
  expect_named(res, c("est", "lci", "uci"))
})

test_that("mad_ratio_boot_cpp: est matches (mad(x)/mad(y))^2", {
  res <- mad_ratio_boot_cpp(x, y, R = 500, seed = 1)
  ref <- (mad(x) / mad(y))^2
  expect_equal(unname(res["est"]), ref, tolerance = 1e-10)
})

test_that("mad_ratio_boot_cpp: est > 0", {
  res <- mad_ratio_boot_cpp(x, y, R = 500, seed = 1)
  expect_gt(unname(res["est"]), 0)
})

test_that("mad_ratio_boot_cpp: identical samples give est = 1", {
  res <- mad_ratio_boot_cpp(x, x, R = 500, seed = 1)
  expect_equal(unname(res["est"]), 1, tolerance = 1e-10)
})
