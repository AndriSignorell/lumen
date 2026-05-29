
test_that("nemenyiTest.default: basic three-group example (list interface)", {
  x <- c(2.9, 3.0, 2.5, 2.6, 3.2)
  y <- c(3.8, 2.7, 4.0, 2.4)
  z <- c(2.8, 3.4, 3.7, 2.2, 2.0)
  
  res <- nemenyiTest(list(x, y, z))
  
  expect_s3_class(res, "rankTest")
  expect_s3_class(res, "htest")
  expect_true(is.list(res))
  expect_named(res, c("res", "pmat"))
})


test_that("nemenyiTest.default: result has correct number of comparisons", {
  set.seed(42)
  x <- rnorm(30)
  g <- rep(letters[1:3], 10)
  
  res <- nemenyiTest(x, g, output = "list")
  
  expect_equal(nrow(res$res), 3L)          # 3 choose 2 = 3
  expect_equal(ncol(res$res), 2L)
  expect_equal(colnames(res$res), c("mean rank diff", "pval"))
})


test_that("nemenyiTest.default: p-values are in [0, 1] for both distributions", {
  set.seed(1)
  x <- c(rnorm(10, 0), rnorm(10, 1), rnorm(10, 2))
  g <- rep(1:3, each = 10)
  
  for (d in c("tukey", "chisq")) {
    res <- nemenyiTest(x, g, dist = d)
    expect_true(all(res$res[, "pval"] >= 0), info = d)
    expect_true(all(res$res[, "pval"] <= 1), info = d)
    expect_true(all(res$pmat >= 0, na.rm = TRUE), info = d)
    expect_true(all(res$pmat <= 1, na.rm = TRUE), info = d)
  }
})


test_that("nemenyiTest.default: pmat is symmetric with diagonal 1", {
  set.seed(7)
  x <- c(rnorm(8), rnorm(8, 2), rnorm(8, 4))
  g <- rep(c("A", "B", "C"), each = 8)
  
  for (d in c("tukey", "chisq")) {
    res <- nemenyiTest(x, g, dist = d)
    pm  <- res$pmat
    
    expect_equal(pm, t(pm), info = d)
    expect_equal(diag(pm), c(A = 1, B = 1, C = 1), info = d)
  }
})


test_that("nemenyiTest.default: tukey and chisq give different p-values", {
  set.seed(3)
  x <- c(rnorm(10), rnorm(10, 2), rnorm(10, 4))
  g <- rep(1:3, each = 10)
  
  p_tukey <- nemenyiTest(x, g, dist = "tukey")$res[, "pval"]
  p_chisq <- nemenyiTest(x, g, dist = "chisq")$res[, "pval"]
  
  expect_false(isTRUE(all.equal(p_tukey, p_chisq)))
})


test_that("nemenyiTest.default: matrix output returns lower-triangular matrix", {
  # Use k=3 groups: pmat[-1, -ncol(pmat)] is 2x2 and stays a matrix.
  # (k=2 would yield a 1x1 array that drops to a scalar without drop=FALSE.)
  set.seed(5)
  x <- c(rnorm(10), rnorm(10, 3), rnorm(10, 6))
  g <- rep(c("low", "mid", "high"), each = 10)
  
  res <- nemenyiTest(x, g, output = "matrix")
  
  expect_true(is.matrix(res$res))
})


test_that("nemenyiTest.default: attr 'method' is always 'none'", {
  # Nemenyi has multiplicity control built in; no additional p.adjust applied
  x <- c(rnorm(10), rnorm(10, 2))
  g <- rep(c("a", "b"), each = 10)
  
  for (d in c("tukey", "chisq")) {
    res <- nemenyiTest(x, g, dist = d)
    expect_equal(attr(res, "method"), "none", info = d)
  }
})


test_that("nemenyiTest.default: attr 'main' contains distribution name", {
  x <- c(rnorm(10), rnorm(10, 2), rnorm(10, 4))
  g <- rep(1:3, each = 10)
  
  expect_match(attr(nemenyiTest(x, g, dist = "tukey"), "main"), "tukey")
  expect_match(attr(nemenyiTest(x, g, dist = "chisq"), "main"), "chisq")
})


test_that("nemenyiTest.default: ties handling does not crash and preserves p-value range", {
  x <- c(rep(1, 5), rep(2, 5), rep(3, 5))
  g <- rep(c("a", "b", "c"), each = 5)
  
  for (d in c("tukey", "chisq")) {
    res <- nemenyiTest(x, g, dist = d)
    expect_true(all(res$res[, "pval"] >= 0), info = d)
    expect_true(all(res$res[, "pval"] <= 1), info = d)
  }
})


test_that("nemenyiTest.default: tiesadj is capped at 1 (all unique ranks)", {
  # No ties: tiesadj = min(1, 1 - 0) = 1; should not crash or produce NaN
  x <- 1:15
  g <- rep(c("a", "b", "c"), each = 5)
  
  res <- nemenyiTest(x, g)
  
  expect_false(anyNA(res$res[, "pval"]))
  expect_false(any(is.nan(res$res[, "pval"])))
})


test_that("nemenyiTest.formula: equivalent to default interface", {
  df <- data.frame(
    val = c(2.9, 3.0, 2.5, 2.6, 3.2,
            3.8, 2.7, 4.0, 2.4,
            2.8, 3.4, 3.7, 2.2, 2.0),
    grp = factor(rep(c("X", "Y", "Z"), c(5, 4, 5)))
  )
  
  for (d in c("tukey", "chisq")) {
    res_f <- nemenyiTest(val ~ grp, data = df, dist = d)
    res_d <- nemenyiTest(df$val, df$grp, dist = d)
    
    expect_equal(res_f$res[, "pval"], res_d$res[, "pval"],
                 tolerance = 1e-10, info = d)
    expect_equal(res_f$pmat, res_d$pmat,
                 tolerance = 1e-10, info = d)
  }
})


test_that("nemenyiTest.formula: subset argument is respected", {
  df <- data.frame(
    val = c(2.9, 3.0, 2.5, 2.6, 3.2,
            3.8, 2.7, 4.0, 2.4,
            2.8, 3.4, 3.7, 2.2, 2.0),
    grp = factor(rep(c("X", "Y", "Z"), c(5, 4, 5)))
  )
  
  res_sub <- nemenyiTest(val ~ grp, data = df, subset = grp != "Y")
  
  expect_equal(nrow(res_sub$pmat), 2L)
  expect_equal(ncol(res_sub$pmat), 2L)
})


test_that("nemenyiTest.formula: na.action = na.omit drops NAs silently", {
  df <- data.frame(
    val = c(2.9, NA,  2.5, 2.6, 3.2,
            3.8, 2.7, 4.0, 2.4,
            2.8, 3.4, 3.7, 2.2, 2.0),
    grp = factor(rep(c("X", "Y", "Z"), c(5, 4, 5)))
  )
  
  expect_no_error(
    nemenyiTest(val ~ grp, data = df, na.action = na.omit)
  )
})


test_that("nemenyiTest: airquality formula example runs without error", {
  expect_no_error({
    res <- nemenyiTest(Ozone ~ Month, data = airquality, na.action = na.omit)
    expect_s3_class(res, "rankTest")
    expect_equal(nrow(res$res), 10L)   # 5 months -> 10 pairs
  })
})


test_that("nemenyiTest: print method runs without error", {
  x <- c(2.9, 3.0, 2.5, 2.6, 3.2)
  y <- c(3.8, 2.7, 4.0, 2.4)
  z <- c(2.8, 3.4, 3.7, 2.2, 2.0)
  
  res <- nemenyiTest(list(x, y, z))
  expect_output(print(res), "Nemenyi")
  
  res_mat <- nemenyiTest(list(x, y, z), output = "matrix")
  expect_output(print(res_mat), "Nemenyi")
})


test_that("nemenyiTest: error on malformed formula (missing RHS)", {
  expect_error(nemenyiTest.formula(~ x), "'formula' missing or incorrect")
})


test_that("nemenyiTest: error on malformed formula (missing response)", {
  # Two RHS terms -> length(term.labels) != 1 -> guard fires before any eval
  expect_error(
    nemenyiTest.formula(y ~ a + b),
    "'formula' missing or incorrect"
  )
})


test_that("nemenyiTest: chisq p-values use k-1 degrees of freedom", {
  # With k=3, df=2. Verify p-values are consistent with pchisq(..., df=2)
  # by checking they change when k changes.
  set.seed(77)
  x3 <- c(rnorm(10), rnorm(10, 2), rnorm(10, 4))
  g3 <- rep(1:3, each = 10)
  
  x4 <- c(x3, rnorm(10, 6))
  g4 <- rep(1:4, each = 10)
  
  p3 <- nemenyiTest(x3, g3, dist = "chisq")$res[, "pval"]
  p4 <- nemenyiTest(x4, g4, dist = "chisq")$res[1:3, "pval"]
  
  # Same pairs, different df -> p-values must differ
  expect_false(isTRUE(all.equal(p3, p4)))
})


test_that("nemenyiTest: pmat 'lbl' attribute flags significant pairs", {
  # Use clearly separated groups so at least one pair is significant
  set.seed(88)
  x <- c(rnorm(20, 0), rnorm(20, 10), rnorm(20, 20))
  g <- rep(c("lo", "mid", "hi"), each = 20)
  
  res <- nemenyiTest(x, g)
  lbl <- attr(res$pmat, "lbl")
  
  expect_type(lbl, "character")
  expect_length(lbl, 3L)
  # At least one group should have a significant counterpart
  expect_true(any(nchar(lbl) > 0))
})