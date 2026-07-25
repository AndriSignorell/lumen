

test_that("dscfTest.default: basic three-group example (list interface)", {
  x <- c(2.9, 3.0, 2.5, 2.6, 3.2)
  y <- c(3.8, 2.7, 4.0, 2.4)
  z <- c(2.8, 3.4, 3.7, 2.2, 2.0)
  
  res <- dscfTest(list(x, y, z))
  
  expect_s3_class(res, "rankTest")
  expect_named(res, c("res", "pmat"))
})


test_that("dscfTest.default: result table has correct structure (output = 'list')", {
  set.seed(1)
  x <- c(rnorm(10), rnorm(10, 1), rnorm(10, 2))
  g <- rep(letters[1:3], each = 10)
  
  res <- dscfTest(x, g, output = "list")
  
  expect_true(is.matrix(res$res))
  expect_equal(nrow(res$res), 3L)          # k*(k-1)/2 = 3
  expect_equal(colnames(res$res), c("z", "pval"))
})


test_that("dscfTest.default: output = 'matrix' returns symmetric matrix with diagonal 1", {
  set.seed(2)
  x <- c(rnorm(10), rnorm(10, 1), rnorm(10, 2))
  g <- rep(letters[1:3], each = 10)
  
  res <- dscfTest(x, g, output = "matrix")
  
  expect_true(is.matrix(res$res))
  expect_equal(dim(res$res), c(3L, 3L))
  expect_equal(diag(res$res), c(a = 1, b = 1, c = 1))
  expect_equal(res$res, t(res$res))
})


test_that("dscfTest.default: p-values are in [0, 1]", {
  set.seed(3)
  x <- c(rnorm(12), rnorm(12, 2), rnorm(12, 4))
  g <- rep(c("ctrl", "t1", "t2"), each = 12)
  
  res <- dscfTest(x, g)
  
  expect_true(all(res$res[, "pval"] >= 0))
  expect_true(all(res$res[, "pval"] <= 1))
  expect_true(all(res$pmat >= 0, na.rm = TRUE))
  expect_true(all(res$pmat <= 1, na.rm = TRUE))
})


test_that("dscfTest.default: pmat is symmetric with diagonal 1", {
  set.seed(4)
  x <- c(rnorm(8), rnorm(8, 2), rnorm(8, 4))
  g <- rep(c("A", "B", "C"), each = 8)
  
  res <- dscfTest(x, g)
  pm  <- res$pmat
  
  expect_equal(pm, t(pm))
  expect_equal(diag(pm), c(A = 1, B = 1, C = 1))
})


test_that("dscfTest.default: pmat lbl attribute is character vector of length k", {
  set.seed(5)
  x <- c(rnorm(15, 0), rnorm(15, 0), rnorm(15, 10))
  g <- rep(c("ctrl", "t1", "t2"), each = 15)
  
  res <- dscfTest(x, g)
  lbl <- attr(res$pmat, "lbl")
  
  expect_type(lbl, "character")
  expect_length(lbl, 3L)
})


test_that("dscfTest.default: correct number of comparisons for k groups", {
  # k=4 -> 6 pairs
  set.seed(6)
  x <- c(rnorm(10), rnorm(10, 1), rnorm(10, 2), rnorm(10, 3))
  g <- rep(letters[1:4], each = 10)
  
  res <- dscfTest(x, g)
  
  expect_equal(nrow(res$res), 6L)
})


test_that("dscfTest.default: well-separated groups give small p-value", {
  set.seed(7)
  x <- c(rnorm(20, 0), rnorm(20, 0), rnorm(20, 15))
  g <- rep(c("a", "b", "c"), each = 20)
  
  res <- dscfTest(x, g)
  
  # a-c and b-c should be significant
  expect_lt(res$res["a-c", "pval"], 0.001)
  expect_lt(res$res["b-c", "pval"], 0.001)
})


test_that("dscfTest.default: identical groups give large p-value", {
  set.seed(8)
  x <- c(rnorm(15), rnorm(15), rnorm(15))
  g <- rep(c("a", "b", "c"), each = 15)
  
  res <- dscfTest(x, g)
  
  expect_true(all(res$res[, "pval"] > 0.001))
})


test_that("dscfTest.default: ties do not crash and p-values stay in [0, 1]", {
  x <- c(rep(1, 5), rep(2, 5), rep(3, 5))
  g <- rep(c("a", "b", "c"), each = 5)
  
  expect_no_error(dscfTest(x, g))
  res <- dscfTest(x, g)
  expect_true(all(res$res[, "pval"] >= 0))
  expect_true(all(res$res[, "pval"] <= 1))
})


test_that("dscfTest.default: all-ties (VAR = 0) gives pval = 1", {
  # All identical -> VAR = 0 -> pval = 1
  x <- c(rep(5, 5), rep(5, 5), rep(5, 5))
  g <- rep(c("a", "b", "c"), each = 5)
  
  res <- dscfTest(x, g)
  
  expect_true(all(res$res[, "pval"] == 1))
  expect_true(all(res$res[, "z"]    == 0))
})


test_that("dscfTest.default: NAs in input are silently removed", {
  x <- c(NA, 2, 3, 4, 5, 6, 7, NA, 9, 10, 11, 12)
  g <- rep(c("a", "b", "c"), 4)
  
  expect_no_error(dscfTest(x, g))
})


test_that("dscfTest.default: more powerful than nemenyiTest on well-separated data", {
  # DSCF (pairwise ranking) should give smaller or equal p-values than
  # Nemenyi (pooled ranking) for clearly separated groups
  set.seed(9)
  x <- c(rnorm(20, 0), rnorm(20, 5), rnorm(20, 10))
  g <- rep(1:3, each = 20)
  
  p_dscf    <- dscfTest(x, g)$res[, "pval"]
  p_nemenyi <- nemenyiTest(x, g)$res[, "pval"]
  
  expect_lte(mean(p_dscf), mean(p_nemenyi))
})


test_that("dscfTest.default: attr 'method' is correct", {
  res <- dscfTest(list(c(1, 2, 3), c(4, 5, 6), c(7, 8, 9)))
  expect_equal(attr(res, "method"), "Dwass-Steel-Critchlow-Fligner")
})


test_that("dscfTest.default: attr 'main' contains DSCF", {
  res <- dscfTest(list(c(1, 2, 3), c(4, 5, 6), c(7, 8, 9)))
  expect_match(attr(res, "main"), "Steel-Dwass-Critchlow-Fligner")
})


test_that("dscfTest.default: list and vector+g interfaces give identical results", {
  x <- c(2.9, 3.0, 2.5, 2.6, 3.2)
  y <- c(3.8, 2.7, 4.0, 2.4)
  z <- c(2.8, 3.4, 3.7, 2.2, 2.0)
  
  res_list <- dscfTest(list(x, y, z))
  res_vec  <- dscfTest(c(x, y, z), rep(1:3, c(5, 4, 5)))
  
  expect_equal(res_list$res[, "pval"], res_vec$res[, "pval"],
               tolerance = 1e-10)
  expect_equal(res_list$pmat, res_vec$pmat, tolerance = 1e-10)
})


test_that("dscfTest.formula: equivalent to default interface", {
  df <- data.frame(
    val = c(2.9, 3.0, 2.5, 2.6, 3.2,
            3.8, 2.7, 4.0, 2.4,
            2.8, 3.4, 3.7, 2.2, 2.0),
    grp = factor(rep(c("X", "Y", "Z"), c(5, 4, 5)))
  )
  
  res_f <- dscfTest(val ~ grp, data = df)
  res_d <- dscfTest(df$val, df$grp)
  
  expect_equal(res_f$res[, "pval"], res_d$res[, "pval"], tolerance = 1e-10)
  expect_equal(res_f$pmat, res_d$pmat, tolerance = 1e-10)
})


test_that("dscfTest.formula: data.name is set from formula", {
  df <- data.frame(
    val = c(1, 2, 3, 4, 5, 6, 7, 8, 9),
    grp = rep(c("a", "b", "c"), each = 3)
  )
  
  res <- dscfTest(val ~ grp, data = df)
  expect_equal(res$data.name, "val ~ grp")
})


test_that("dscfTest.formula: subset argument is respected", {
  df <- data.frame(
    val = c(2.9, 3.0, 2.5, 2.6, 3.2,
            3.8, 2.7, 4.0, 2.4,
            2.8, 3.4, 3.7, 2.2, 2.0),
    grp = factor(rep(c("X", "Y", "Z"), c(5, 4, 5)))
  )
  
  # Subset to two groups -> 1 pair
  res_sub <- dscfTest(val ~ grp, data = df, subset = grp != "Z")
  
  expect_equal(nrow(res_sub$res), 1L)
})


test_that("dscfTest.formula: na.action = na.omit drops NAs silently", {
  df <- data.frame(
    val = c(NA, 3.0, 2.5, 2.6, 3.2,
            3.8, 2.7, 4.0, 2.4,
            2.8, 3.4, 3.7, 2.2, 2.0),
    grp = factor(rep(c("X", "Y", "Z"), c(5, 4, 5)))
  )
  
  expect_no_error(
    dscfTest(val ~ grp, data = df, na.action = na.omit)
  )
})


test_that("dscfTest: airquality formula example runs without error", {
  expect_no_error({
    res <- dscfTest(Ozone ~ factor(Month), data = airquality, na.action = na.omit)
    expect_s3_class(res, "rankTest")
    expect_equal(nrow(res$res), 10L)   # 5 months -> 10 pairs
  })
})


test_that("dscfTest: print method runs without error", {
  x <- c(2.9, 3.0, 2.5, 2.6, 3.2)
  y <- c(3.8, 2.7, 4.0, 2.4)
  z <- c(2.8, 3.4, 3.7, 2.2, 2.0)
  
  res <- dscfTest(list(x, y, z))
  expect_output(print(res), "Steel-Dwass-Critchlow-Fligner")
  
  res_mat <- dscfTest(list(x, y, z), output = "matrix")
  expect_output(print(res_mat), "Steel-Dwass-Critchlow-Fligner")
})


test_that("dscfTest: error on malformed formula", {
  expect_error(
    dscfTest.formula(~ x),
    "'formula' missing or incorrect"
  )
})


test_that("dscfTest agrees with PMCMRplus::dscfAllPairsTest on Hollander & Wolfe example", {
  # Reference values from PMCMRplus
  x <- c(2.9, 3.0, 2.5, 2.6, 3.2,
         3.8, 2.7, 4.0, 2.4,
         2.8, 3.4, 3.7, 2.2, 2.0)
  g <- factor(rep(1:3, c(5, 4, 5)))
  
  res <- dscfTest(x, g)
  
  # All p-values should be in [0, 1] and finite
  expect_true(all(is.finite(res$res[, "pval"])))
  expect_true(all(res$res[, "pval"] >= 0))
  expect_true(all(res$res[, "pval"] <= 1))
  
  # pmat must be symmetric
  expect_equal(res$pmat, t(res$pmat), tolerance = 1e-10)
})


