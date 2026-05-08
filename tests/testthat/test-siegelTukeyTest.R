

# tests/testthat/test-siegelTukeyTest.R

test_that("siegelTukeyTest.default returns htest object", {
  x <- c(1, 2, 3, 4, 5)
  y <- c(1, 3, 5, 7, 9)
  
  res <- siegelTukeyTest(x, y)
  
  expect_s3_class(res, "htest")
  expect_named(res, c(
    "statistic", "parameter", "p.value", "null.value",
    "alternative", "method", "data.name", "exact", "ties"
  ))
  expect_named(res$statistic, "W")
  expect_true(is.numeric(res$statistic))
  expect_true(is.numeric(res$p.value))
  expect_gte(res$p.value, 0)
  expect_lte(res$p.value, 1)
  expect_equal(res$method, "Siegel-Tukey test for scale differences")
})


test_that("siegelTukeyTest requires y", {
  x <- 1:5
  
  expect_error(
    siegelTukeyTest(x),
    "'y' is missing",
    fixed = TRUE
  )
})


test_that("alternative argument is validated", {
  x <- 1:5
  y <- 6:10
  
  expect_error(
    siegelTukeyTest(x, y, alternative = "invalid")
  )
  
  expect_equal(
    siegelTukeyTest(x, y, alternative = "less")$alternative,
    "less"
  )
  
  expect_equal(
    siegelTukeyTest(x, y, alternative = "greater")$alternative,
    "greater"
  )
})


test_that("mu must be a single finite number", {
  x <- 1:5
  y <- 6:10
  
  expect_error(
    siegelTukeyTest(x, y, mu = c(0, 1)),
    "'mu' must be a single number",
    fixed = TRUE
  )
  
  expect_error(
    siegelTukeyTest(x, y, mu = Inf),
    "'mu' must be a single number",
    fixed = TRUE
  )
  
  expect_equal(
    siegelTukeyTest(x, y, mu = 2)$null.value,
    c(mu = 2)
  )
})


test_that("formula interface works for two independent groups", {
  dat <- data.frame(
    value = c(1, 2, 3, 4, 1, 3, 5, 7),
    group = rep(c("a", "b"), each = 4)
  )
  
  res_formula <- siegelTukeyTest(value ~ group, data = dat)
  res_default <- siegelTukeyTest(
    x = dat$value[dat$group == "a"],
    y = dat$value[dat$group == "b"]
  )
  
  expect_s3_class(res_formula, "htest")
  expect_equal(res_formula$p.value, res_default$p.value)
  expect_equal(res_formula$statistic, res_default$statistic)
})


test_that("formula interface rejects incorrect formula", {
  dat <- data.frame(
    value = 1:5,
    group = rep(c("a", "b"), length.out = 5)
  )
  
  expect_error(
    siegelTukeyTest(~ value, data = dat),
    "'formula' missing or incorrect",
    fixed = TRUE
  )
})


test_that("adjust.median changes result when medians differ", {
  x <- c(1, 2, 3, 4, 5)
  y <- c(10, 12, 14, 16, 18)
  
  res_unadjusted <- siegelTukeyTest(x, y, adjust.median = FALSE)
  res_adjusted   <- siegelTukeyTest(x, y, adjust.median = TRUE)
  
  expect_s3_class(res_adjusted, "htest")
  expect_false(isTRUE(all.equal(
    res_unadjusted$statistic,
    res_adjusted$statistic
  )))
})


test_that("ties are detected", {
  x <- c(1, 2, 2, 3, 4)
  y <- c(1, 2, 3, 3, 5)
  
  res <- siegelTukeyTest(x, y)
  
  expect_true(res$ties)
  expect_false(res$exact)
})


test_that("exact p-value is used for small samples without ties", {
  x <- c(1, 2, 3, 4)
  y <- c(5, 6, 7, 8)
  
  res <- siegelTukeyTest(x, y)
  
  expect_true(res$exact)
  expect_false(res$ties)
  expect_gte(res$p.value, 0)
  expect_lte(res$p.value, 1)
})


test_that("normal approximation is used for large samples", {
  x <- 1:60
  y <- 101:160
  
  res <- siegelTukeyTest(x, y)
  
  expect_false(res$exact)
  expect_false(res$ties)
  expect_gte(res$p.value, 0)
  expect_lte(res$p.value, 1)
})


test_that(".siegelTukeyRank returns expected columns and length for even n", {
  x <- c(1, 2, 3, 4, 5, 6)
  g <- c(0, 0, 0, 1, 1, 1)
  
  res <- .siegelTukeyRank(x, g)
  
  expect_s3_class(res, "data.frame")
  expect_named(res, c("sort.x", "sort.id", "unique.ranks", "raw.ranks"))
  expect_equal(nrow(res), length(x))
  expect_equal(res$sort.x, sort(x))
})


test_that(".siegelTukeyRank drops one median observation for odd n", {
  x <- c(1, 2, 3, 4, 5)
  g <- c(0, 0, 1, 1, 1)
  
  res <- .siegelTukeyRank(x, g, dropMedian = TRUE)
  
  expect_equal(nrow(res), length(x) - 1)
  expect_false(3 %in% res$sort.x)
})


test_that(".siegelTukeyRank averages ranks for tied x values", {
  x <- c(1, 2, 2, 3)
  g <- c(0, 0, 1, 1)
  
  res <- .siegelTukeyRank(x, g, dropMedian = FALSE)
  
  tied_ranks <- res$unique.ranks[res$sort.x == 2]
  
  expect_equal(length(unique(tied_ranks)), 1)
})