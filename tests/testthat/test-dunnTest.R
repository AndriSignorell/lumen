
test_that("dunnTest.default: basic three-group example (list interface)", {
  x <- c(2.9, 3.0, 2.5, 2.6, 3.2)
  y <- c(3.8, 2.7, 4.0, 2.4)
  z <- c(2.8, 3.4, 3.7, 2.2, 2.0)
  
  res <- dunnTest(list(x, y, z))
  
  expect_s3_class(res, "rankTest")
  expect_s3_class(res, "htest")
  expect_true(is.list(res))
  expect_named(res, c("res", "pmat"))
})


test_that("dunnTest.default: result has correct number of comparisons", {
  # k groups -> k*(k-1)/2 pairwise comparisons
  set.seed(42)
  x <- rnorm(30)
  g <- rep(letters[1:3], 10)
  
  res <- dunnTest(x, g, output = "list")
  
  expect_equal(nrow(res$res), 3L)          # 3 choose 2 = 3
  expect_equal(ncol(res$res), 2L)          # mean rank diff + pval
  expect_equal(colnames(res$res), c("mean rank diff", "pval"))
})


test_that("dunnTest.default: p-values are in [0, 1]", {
  set.seed(1)
  x <- c(rnorm(10, 0), rnorm(10, 1), rnorm(10, 2))
  g <- rep(1:3, each = 10)
  
  res <- dunnTest(x, g)
  
  expect_true(all(res$res[, "pval"] >= 0))
  expect_true(all(res$res[, "pval"] <= 1))
  expect_true(all(res$pmat >= 0, na.rm = TRUE))
  expect_true(all(res$pmat <= 1, na.rm = TRUE))
})


test_that("dunnTest.default: pmat is symmetric with diagonal 1", {
  set.seed(7)
  x <- c(rnorm(8), rnorm(8, 2), rnorm(8, 4))
  g <- rep(c("A", "B", "C"), each = 8)
  
  res <- dunnTest(x, g)
  pm  <- res$pmat
  
  expect_equal(pm, t(pm))                         # symmetric
  expect_equal(diag(pm), c(A = 1, B = 1, C = 1)) # diagonal = 1
})


test_that("dunnTest.default: matrix output returns lower-triangular matrix", {
  set.seed(3)
  x <- c(rnorm(10), rnorm(10, 3))
  g <- rep(c("low", "high"), each = 10)
  
  res <- dunnTest(x, g, output = "matrix")
  
  # For k=2 the result should be a 1×1 matrix (lower tri without last col)
  expect_true(is.matrix(res$res))
})


test_that("dunnTest.default: p-adjustment methods produce different results", {
  set.seed(99)
  x <- c(rnorm(12, 0), rnorm(12, 1), rnorm(12, 2), rnorm(12, 3))
  g <- rep(1:4, each = 12)
  
  p_none <- dunnTest(x, g, method = "none")$res[, "pval"]
  p_bonf <- dunnTest(x, g, method = "bonferroni")$res[, "pval"]
  p_holm <- dunnTest(x, g, method = "holm")$res[, "pval"]
  
  # Bonferroni and Holm must be >= unadjusted
  expect_true(all(p_bonf >= p_none - .Machine$double.eps))
  expect_true(all(p_holm >= p_none - .Machine$double.eps))
  
  # Bonferroni and Holm are not identical in general (>= 4 groups)
  expect_false(identical(p_bonf, p_holm))
})


test_that("dunnTest.default: alternative = 'less' / 'greater' give one-sided p-values", {
  set.seed(5)
  x <- c(rnorm(10, 0), rnorm(10, 3))
  g <- rep(c("low", "high"), each = 10)
  
  p_two   <- dunnTest(x, g, method = "none", alternative = "two.sided")$res[, "pval"]
  p_less  <- dunnTest(x, g, method = "none", alternative = "less")$res[, "pval"]
  p_great <- dunnTest(x, g, method = "none", alternative = "greater")$res[, "pval"]
  
  # One-sided p-values should differ from two-sided
  expect_false(isTRUE(all.equal(p_two, p_less)))
  expect_false(isTRUE(all.equal(p_two, p_great)))
  
  # less + greater should roughly sum to 1 (before adjustment, single comparison)
  expect_equal(p_less + p_great, rep(1, length(p_less)), tolerance = 1e-10)
})


test_that("dunnTest.formula: equivalent to default interface", {
  df <- data.frame(
    val = c(2.9, 3.0, 2.5, 2.6, 3.2,
            3.8, 2.7, 4.0, 2.4,
            2.8, 3.4, 3.7, 2.2, 2.0),
    grp = factor(rep(c("X", "Y", "Z"), c(5, 4, 5)))
  )
  
  res_formula <- dunnTest(val ~ grp, data = df)
  res_default <- dunnTest(df$val, df$grp)
  
  expect_equal(res_formula$res[, "pval"], res_default$res[, "pval"],
               tolerance = 1e-10)
  expect_equal(res_formula$pmat, res_default$pmat,
               tolerance = 1e-10)
})


test_that("dunnTest.formula: subset argument is respected", {
  df <- data.frame(
    val = c(2.9, 3.0, 2.5, 2.6, 3.2,
            3.8, 2.7, 4.0, 2.4,
            2.8, 3.4, 3.7, 2.2, 2.0),
    grp = factor(rep(c("X", "Y", "Z"), c(5, 4, 5)))
  )
  
  # Only groups X and Z
  res_sub <- dunnTest(val ~ grp, data = df, subset = grp != "Y")
  
  expect_equal(nrow(res_sub$pmat), 2L)
  expect_equal(ncol(res_sub$pmat), 2L)
})


test_that("dunnTest.formula: na.action = na.omit drops NAs silently", {
  df <- data.frame(
    val = c(2.9, NA,  2.5, 2.6, 3.2,
            3.8, 2.7, 4.0, 2.4,
            2.8, 3.4, 3.7, 2.2, 2.0),
    grp = factor(rep(c("X", "Y", "Z"), c(5, 4, 5)))
  )
  
  expect_no_error(
    dunnTest(val ~ grp, data = df, na.action = na.omit)
  )
})


test_that("dunnTest: airquality formula example runs without error", {
  expect_no_error({
    res <- dunnTest(Ozone ~ Month, data = airquality, na.action = na.omit)
    expect_s3_class(res, "rankTest")
    # 5 months -> 10 pairwise comparisons
    expect_equal(nrow(res$res), 10L)
  })
})


test_that("dunnTest: print method runs without error", {
  x <- c(2.9, 3.0, 2.5, 2.6, 3.2)
  y <- c(3.8, 2.7, 4.0, 2.4)
  z <- c(2.8, 3.4, 3.7, 2.2, 2.0)
  
  res <- dunnTest(list(x, y, z))
  expect_output(print(res), "Dunn")
  
  res_mat <- dunnTest(list(x, y, z), output = "matrix")
  expect_output(print(res_mat), "Dunn")
})


test_that("dunnTest: error on malformed formula", {
  expect_error(dunnTest.formula(~ x), "'formula' missing or incorrect")
})


test_that("dunnTest: ties handling does not crash and preserves p-value range", {
  # Lots of ties
  x <- c(rep(1, 5), rep(2, 5), rep(3, 5))
  g <- rep(c("a", "b", "c"), each = 5)
  
  res <- dunnTest(x, g, method = "none")
  
  expect_true(all(res$res[, "pval"] >= 0))
  expect_true(all(res$res[, "pval"] <= 1))
})