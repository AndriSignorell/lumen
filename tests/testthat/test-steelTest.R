
test_that("steelTest.default: basic three-group example (list interface)", {
  x <- c(2.9, 3.0, 2.5, 2.6, 3.2)
  y <- c(3.8, 2.7, 4.0, 2.4)
  z <- c(2.8, 3.4, 3.7, 2.2, 2.0)
  
  res <- steelTest(list(x, y, z))
  
  expect_s3_class(res, "rankTest")
  expect_s3_class(res, "htest")
  expect_named(res, c("statistic", "p.value", "res", "pmat", "corr"))
})


test_that("steelTest.default: result table has correct structure (output = 'list')", {
  set.seed(1)
  x <- c(rnorm(10), rnorm(10, 1), rnorm(10, 2))
  g <- rep(c("ctrl", "t1", "t2"), each = 10)
  
  res <- steelTest(x, g, control = "ctrl", output = "list")
  
  expect_true(is.matrix(res$res))
  expect_equal(nrow(res$res), 2L)          # m = k - 1 = 2 treatments
  expect_equal(colnames(res$res), c("W", "z", "pval"))
})


test_that("steelTest.default: output = 'matrix' returns a full group x group matrix", {
  set.seed(2)
  x <- c(rnorm(10), rnorm(10, 1), rnorm(10, 2))
  g <- rep(c("ctrl", "t1", "t2"), each = 10)
  
  res <- steelTest(x, g, control = "ctrl", output = "matrix")
  
  expect_true(is.matrix(res$res))
  expect_equal(dim(res$res), c(3L, 3L))
})


test_that("steelTest.default: p-values are in [0, 1]", {
  set.seed(3)
  x <- c(rnorm(12), rnorm(12, 2), rnorm(12, 4))
  g <- rep(c("ctrl", "t1", "t2"), each = 12)
  
  res <- steelTest(x, g, control = "ctrl")
  
  expect_true(all(res$res[, "pval"] >= 0))
  expect_true(all(res$res[, "pval"] <= 1))
  expect_gte(res$p.value, 0)
  expect_lte(res$p.value, 1)
})


test_that("steelTest.default: global p-value <= min pairwise p-value", {
  # Global statistic is max(|z|); its p-value should be <= any single pairwise
  set.seed(4)
  x <- c(rnorm(15), rnorm(15, 3), rnorm(15, 6))
  g <- rep(c("ctrl", "t1", "t2"), each = 15)
  
  res <- steelTest(x, g, control = "ctrl")
  
  expect_lte(res$p.value, min(res$res[, "pval"]) + 1e-10)
})


test_that("steelTest.default: pmat diagonal is 1, off-diagonal NA except treatment-control", {
  set.seed(5)
  x <- c(rnorm(10), rnorm(10, 1), rnorm(10, 2))
  g <- rep(c("ctrl", "t1", "t2"), each = 10)
  
  res <- steelTest(x, g, control = "ctrl")
  pm  <- res$pmat
  
  expect_equal(diag(pm), c(ctrl = 1, t1 = 1, t2 = 1))
  expect_true(all(!is.na(pm["t1", "ctrl"])))
  expect_true(all(!is.na(pm["t2", "ctrl"])))
  expect_true(is.na(pm["ctrl", "t1"]))
  expect_true(is.na(pm["t1", "t2"]))
})


test_that("steelTest.default: pmat lbl attribute is a character vector of length k", {
  set.seed(6)
  x <- c(rnorm(15, 0), rnorm(15, 0), rnorm(15, 10))
  g <- rep(c("ctrl", "t1", "t2"), each = 15)
  
  res <- steelTest(x, g, control = "ctrl")
  lbl <- attr(res$pmat, "lbl")
  
  expect_type(lbl, "character")
  expect_length(lbl, 3L)
})


test_that("steelTest.default: well-separated treatment gives small p-value", {
  set.seed(7)
  x <- c(rnorm(20, 0), rnorm(20, 0), rnorm(20, 15))
  g <- rep(c("ctrl", "t1", "t2"), each = 20)
  
  res <- steelTest(x, g, control = "ctrl")
  
  expect_lt(res$res["t2-ctrl", "pval"], 0.001)
  expect_lt(res$p.value, 0.001)
})


test_that("steelTest.default: identical groups give large p-value", {
  set.seed(8)
  x <- c(rnorm(15), rnorm(15), rnorm(15))
  g <- rep(c("ctrl", "t1", "t2"), each = 15)
  
  res <- steelTest(x, g, control = "ctrl")
  
  expect_gt(res$p.value, 0.001)
})


test_that("steelTest.default: default control = first factor level", {
  set.seed(9)
  x <- c(rnorm(10), rnorm(10, 1), rnorm(10, 2))
  g <- factor(rep(c("a", "b", "c"), each = 10))
  
  res_default  <- steelTest(x, g)
  res_explicit <- steelTest(x, g, control = "a")
  
  expect_equal(res_default$res, res_explicit$res, tolerance = 1e-10)
  expect_equal(res_default$p.value, res_explicit$p.value, tolerance = 1e-10)
})


test_that("steelTest.default: invalid control throws informative error", {
  x <- c(rnorm(10), rnorm(10))
  g <- rep(c("a", "b"), each = 10)
  
  expect_error(
    steelTest(x, g, control = "z"),
    "not found in grouping variable"
  )
})


test_that("steelTest.default: m = 1 (single treatment) works", {
  set.seed(10)
  x <- c(rnorm(10), rnorm(10, 2))
  g <- rep(c("ctrl", "trt"), each = 10)
  
  expect_no_error(steelTest(x, g, control = "ctrl"))
  res <- steelTest(x, g, control = "ctrl")
  expect_equal(nrow(res$res), 1L)
  # correlation matrix must be 1x1
  expect_equal(dim(res$corr), c(1L, 1L))
})


test_that("steelTest.default: alternative = 'less' / 'greater' give different p-values", {
  set.seed(11)
  x <- c(rnorm(12, 0), rnorm(12, 3), rnorm(12, 3))
  g <- rep(c("ctrl", "t1", "t2"), each = 12)
  
  p_two   <- steelTest(x, g, control = "ctrl", alternative = "two.sided")$p.value
  p_less  <- steelTest(x, g, control = "ctrl", alternative = "less")$p.value
  p_great <- steelTest(x, g, control = "ctrl", alternative = "greater")$p.value
  
  expect_false(isTRUE(all.equal(p_two, p_less)))
  expect_false(isTRUE(all.equal(p_two, p_great)))
  expect_false(isTRUE(all.equal(p_less, p_great)))
})


test_that("steelTest.default: attr 'alternative' is stored correctly", {
  x <- c(rnorm(10), rnorm(10, 1), rnorm(10, 2))
  g <- rep(c("ctrl", "t1", "t2"), each = 10)
  
  for (alt in c("two.sided", "less", "greater")) {
    res <- steelTest(x, g, control = "ctrl", alternative = alt)
    expect_equal(attr(res, "alternative"), alt, info = alt)
  }
})


test_that("steelTest.default: attr 'method' is 'Steel'", {
  res <- steelTest(list(c(1, 2, 3), c(4, 5, 6), c(7, 8, 9)))
  expect_equal(attr(res, "method"), "Steel")
})


test_that("steelTest.default: correlation matrix is symmetric with unit diagonal", {
  set.seed(12)
  x <- c(rnorm(10), rnorm(10, 1), rnorm(10, 2), rnorm(10, 3))
  g <- rep(c("ctrl", "t1", "t2", "t3"), each = 10)
  
  res <- steelTest(x, g, control = "ctrl")
  R   <- res$corr
  
  expect_equal(R, t(R), tolerance = 1e-10)
  expect_equal(diag(R), rep(1, 3))
  # off-diagonal correlations in (0, 1) for positive group sizes
  expect_true(all(R[lower.tri(R)] > 0))
  expect_true(all(R[lower.tri(R)] < 1))
})


test_that("steelTest.default: ties do not crash and p-values stay in [0, 1]", {
  x <- c(rep(1, 5), rep(2, 5), rep(3, 5))
  g <- rep(c("ctrl", "t1", "t2"), each = 5)
  
  expect_no_error(steelTest(x, g, control = "ctrl"))
  res <- steelTest(x, g, control = "ctrl")
  expect_true(all(res$res[, "pval"] >= 0))
  expect_true(all(res$res[, "pval"] <= 1))
})


test_that("steelTest.default: NAs in input are silently removed", {
  x <- c(NA, 2, 3, 4, 5, 6, 7, NA, 9, 10, 11, 12)
  g <- rep(c("ctrl", "t1", "t2"), 4)
  
  expect_no_error(steelTest(x, g))
})


test_that("steelTest.formula: equivalent to default interface", {
  df <- data.frame(
    val = c(2.9, 3.0, 2.5, 2.6, 3.2,
            3.8, 2.7, 4.0, 2.4,
            2.8, 3.4, 3.7, 2.2, 2.0),
    grp = factor(rep(c("ctrl", "t1", "t2"), c(5, 4, 5)))
  )
  
  res_f <- steelTest(val ~ grp, data = df, control = "ctrl")
  res_d <- steelTest(df$val, df$grp, control = "ctrl")
  
  expect_equal(res_f$res,     res_d$res,     tolerance = 1e-10)
  expect_equal(res_f$p.value, res_d$p.value, tolerance = 1e-10)
})


test_that("steelTest.formula: data.name is set from formula", {
  df <- data.frame(
    val = c(1, 2, 3, 4, 5, 6, 7, 8, 9),
    grp = rep(c("ctrl", "t1", "t2"), each = 3)
  )
  
  res <- steelTest(val ~ grp, data = df)
  expect_equal(res$data.name, "val ~ grp")
})


test_that("steelTest.formula: subset argument is respected", {
  df <- data.frame(
    val = c(2.9, 3.0, 2.5, 2.6, 3.2,
            3.8, 2.7, 4.0, 2.4,
            2.8, 3.4, 3.7, 2.2, 2.0),
    grp = factor(rep(c("ctrl", "t1", "t2"), c(5, 4, 5)))
  )
  
  # Subset to ctrl + t1 only -> m = 1
  res_sub <- steelTest(val ~ grp, data = df,
                       subset = grp != "t2", control = "ctrl")
  
  expect_equal(nrow(res_sub$res), 1L)
})


test_that("steelTest.formula: na.action = na.omit drops NAs silently", {
  df <- data.frame(
    val = c(NA, 3.0, 2.5, 2.6, 3.2,
            3.8, 2.7, 4.0, 2.4,
            2.8, 3.4, 3.7, 2.2, 2.0),
    grp = factor(rep(c("ctrl", "t1", "t2"), c(5, 4, 5)))
  )
  
  expect_no_error(
    steelTest(val ~ grp, data = df, na.action = na.omit)
  )
})


test_that("steelTest: airquality example runs without error", {
  expect_no_error({
    res <- steelTest(Ozone ~ Month, data = airquality,
                     na.action = na.omit, control = "5")
    expect_s3_class(res, "rankTest")
    expect_equal(nrow(res$res), 4L)   # 5 months, control = "5" -> 4 rows
  })
})


test_that("steelTest: error on malformed formula", {
  expect_error(
    steelTest.formula(~ x),
    "'formula' missing or incorrect"
  )
})


test_that("steelTest: list and vector+g interfaces give identical results", {
  x <- c(2.9, 3.0, 2.5, 2.6, 3.2)
  y <- c(3.8, 2.7, 4.0, 2.4)
  z <- c(2.8, 3.4, 3.7, 2.2, 2.0)
  
  res_list <- steelTest(list(x, y, z))
  res_vec  <- steelTest(c(x, y, z), rep(1:3, c(5, 4, 5)))
  
  expect_equal(res_list$res[, "pval"], res_vec$res[, "pval"],
               tolerance = 1e-10)
  expect_equal(res_list$p.value, res_vec$p.value, tolerance = 1e-10)
})

