

test_that("dunnettTest.default: basic three-group example (list interface)", {
  x <- c(2.9, 3.0, 2.5, 2.6, 3.2)
  y <- c(3.8, 2.7, 4.0, 2.4)
  z <- c(2.8, 3.4, 3.7, 2.2, 2.0)
  
  res <- dunnettTest(list(x, y, z))
  
  expect_s3_class(res, "PostHocTest")
  expect_true(is.list(res))
  expect_length(res, 1L)   # one control -> one matrix
})


test_that("dunnettTest.default: result matrix has correct columns", {
  set.seed(1)
  x <- c(rnorm(10), rnorm(10, 1), rnorm(10, 2))
  g <- rep(letters[1:3], each = 10)
  
  res <- dunnettTest(x, g)
  
  mat <- res[[1L]]
  expect_true(is.matrix(mat))
  expect_equal(colnames(mat), c("diff", "lwr.ci", "upr.ci", "pval"))
})


test_that("dunnettTest.default: k-1 rows for k groups (default control)", {
  set.seed(2)
  x <- c(rnorm(10), rnorm(10, 1), rnorm(10, 2), rnorm(10, 3))
  g <- rep(letters[1:4], each = 10)
  
  res <- dunnettTest(x, g)
  
  # default control = first level "a"; k-1 = 3 treatment rows
  expect_equal(nrow(res[[1L]]), 3L)
})


test_that("dunnettTest.default: p-values are in [0, 1]", {
  set.seed(3)
  x <- c(rnorm(12), rnorm(12, 2), rnorm(12, 4))
  g <- rep(c("ctrl", "trt1", "trt2"), each = 12)
  
  res <- dunnettTest(x, g, control = "ctrl")
  
  expect_true(all(res[["ctrl"]][, "pval"] >= 0))
  expect_true(all(res[["ctrl"]][, "pval"] <= 1))
})


test_that("dunnettTest.default: CI contains diff when groups are identical", {
  set.seed(4)
  x <- c(rnorm(15, 5), rnorm(15, 5), rnorm(15, 5))
  g <- rep(c("a", "b", "c"), each = 15)
  
  res <- dunnettTest(x, g, control = "a")
  mat <- res[["a"]]
  
  # CI should contain 0 (no real difference)
  expect_true(all(mat[, "lwr.ci"] <= 0))
  expect_true(all(mat[, "upr.ci"] >= 0))
})


test_that("dunnettTest.default: CI is wider for lower conf.level", {
  set.seed(5)
  x <- c(rnorm(10), rnorm(10, 1), rnorm(10, 2))
  g <- rep(c("ctrl", "t1", "t2"), each = 10)
  
  res95 <- dunnettTest(x, g, control = "ctrl", conf.level = 0.95)
  res80 <- dunnettTest(x, g, control = "ctrl", conf.level = 0.80)
  
  width95 <- res95[["ctrl"]][, "upr.ci"] - res95[["ctrl"]][, "lwr.ci"]
  width80 <- res80[["ctrl"]][, "upr.ci"] - res80[["ctrl"]][, "lwr.ci"]
  
  expect_true(all(width95 > width80))
})


test_that("dunnettTest.default: well-separated treatment gives small p-value", {
  set.seed(6)
  x <- c(rnorm(20, 0), rnorm(20, 0), rnorm(20, 10))
  g <- rep(c("ctrl", "trt1", "trt2"), each = 20)
  
  res <- dunnettTest(x, g, control = "ctrl")
  
  # trt2 clearly different from ctrl
  p_trt2 <- res[["ctrl"]]["trt2-ctrl", "pval"]
  expect_lt(p_trt2, 0.001)
})


test_that("dunnettTest.default: explicit control = first level same as default", {
  set.seed(7)
  x <- c(rnorm(10), rnorm(10, 1), rnorm(10, 2))
  g <- factor(rep(c("a", "b", "c"), each = 10))
  
  res_default  <- dunnettTest(x, g)
  res_explicit <- dunnettTest(x, g, control = "a")
  
  expect_equal(res_default[["a"]], res_explicit[["a"]], tolerance = 1e-10)
})


test_that("dunnettTest.default: multiple control levels produce one matrix each", {
  set.seed(8)
  x <- c(rnorm(10), rnorm(10, 1), rnorm(10, 2), rnorm(10, 3))
  g <- rep(letters[1:4], each = 10)
  
  res <- dunnettTest(x, g, control = c("a", "b"))
  
  expect_length(res, 2L)
  expect_named(res, c("a", "b"))
  expect_true(is.matrix(res[["a"]]))
  expect_true(is.matrix(res[["b"]]))
})


test_that("dunnettTest.default: invalid control level throws error", {
  x <- c(rnorm(10), rnorm(10, 1))
  g <- rep(c("a", "b"), each = 10)
  
  expect_error(
    dunnettTest(x, g, control = "z"),
    "not found in grouping variable"
  )
})


test_that("dunnettTest.default: diff equals treatment mean minus control mean", {
  set.seed(9)
  x <- c(rnorm(10, 0), rnorm(10, 3), rnorm(10, 6))
  g <- factor(rep(c("ctrl", "trt1", "trt2"), each = 10))
  
  res  <- dunnettTest(x, g, control = "ctrl")
  mat  <- res[["ctrl"]]
  means <- tapply(x, g, mean)
  
  expected_diff_trt1 <- means["trt1"] - means["ctrl"]
  expected_diff_trt2 <- means["trt2"] - means["ctrl"]
  
  expect_equal(unname(mat["trt1-ctrl", "diff"]),
               unname(expected_diff_trt1), tolerance = 1e-10)
  expect_equal(unname(mat["trt2-ctrl", "diff"]),
               unname(expected_diff_trt2), tolerance = 1e-10)
})


test_that("dunnettTest.default: attr 'method' equals 'Dunnett'", {
  res <- dunnettTest(list(c(1, 2, 3), c(4, 5, 6), c(7, 8, 9)))
  expect_equal(attr(res, "method"), "Dunnett")
})


test_that("dunnettTest.default: attr 'conf.level' is preserved", {
  x <- c(rnorm(10), rnorm(10, 1), rnorm(10, 2))
  g <- rep(c("a", "b", "c"), each = 10)
  
  res <- dunnettTest(x, g, conf.level = 0.90)
  expect_equal(attr(res, "conf.level"), 0.90)
})


test_that("dunnettTest.default: results are reproducible (seed = 5L in qmvt)", {
  set.seed(42)
  x <- c(rnorm(15), rnorm(15, 1), rnorm(15, 2))
  g <- rep(c("ctrl", "trt1", "trt2"), each = 15)
  
  res1 <- dunnettTest(x, g, control = "ctrl")
  res2 <- dunnettTest(x, g, control = "ctrl")
  
  expect_equal(res1[["ctrl"]], res2[["ctrl"]])
})


test_that("dunnettTest.default: NAs in input are silently removed", {
  x <- c(1, NA, 3, 4, 5, 6, 7, 8, NA, 10, 11, 12)
  g <- rep(c("a", "b", "c"), 4)
  
  expect_no_error(dunnettTest(x, g))
})


test_that("dunnettTest.formula: equivalent to default interface", {
  df <- data.frame(
    val = c(2.9, 3.0, 2.5, 2.6, 3.2,
            3.8, 2.7, 4.0, 2.4,
            2.8, 3.4, 3.7, 2.2, 2.0),
    grp = factor(rep(c("X", "Y", "Z"), c(5, 4, 5)))
  )
  
  res_f <- dunnettTest(val ~ grp, data = df, control = "X")
  res_d <- dunnettTest(df$val, df$grp, control = "X")
  
  expect_equal(res_f[["X"]], res_d[["X"]], tolerance = 1e-10)
})


test_that("dunnettTest.formula: data.name is set from formula", {
  df <- data.frame(
    val = c(1, 2, 3, 4, 5, 6, 7, 8, 9),
    grp = rep(c("a", "b", "c"), each = 3)
  )
  
  res <- dunnettTest(val ~ grp, data = df)
  expect_equal(res$data.name, "val ~ grp")
})


test_that("dunnettTest.formula: subset argument is respected", {
  df <- data.frame(
    val = c(2.9, 3.0, 2.5, 2.6, 3.2,
            3.8, 2.7, 4.0, 2.4,
            2.8, 3.4, 3.7, 2.2, 2.0),
    grp = factor(rep(c("X", "Y", "Z"), c(5, 4, 5)))
  )
  
  # Subset to two groups -> one row in result matrix
  res_sub <- dunnettTest(val ~ grp, data = df,
                         subset = grp != "Z", control = "X")
  
  expect_equal(nrow(res_sub[["X"]]), 1L)
})


test_that("dunnettTest.formula: na.action = na.omit drops NAs silently", {
  df <- data.frame(
    val = c(NA, 3.0, 2.5, 2.6, 3.2,
            3.8, 2.7, 4.0, 2.4,
            2.8, 3.4, 3.7, 2.2, 2.0),
    grp = factor(rep(c("X", "Y", "Z"), c(5, 4, 5)))
  )
  
  expect_no_error(
    dunnettTest(val ~ grp, data = df, na.action = na.omit)
  )
})


test_that("dunnettTest: airquality formula example runs without error", {
  expect_no_error({
    res <- dunnettTest(Ozone ~ Month, data = airquality,
                       na.action = na.omit, control = "5")
    expect_s3_class(res, "PostHocTest")
    # 5 months, control = "5" -> 4 treatment rows
    expect_equal(nrow(res[["5"]]), 4L)
  })
})


test_that("dunnettTest: error on malformed formula", {
  expect_error(
    dunnettTest.formula(~ x),
    "'formula' missing or incorrect"
  )
})


test_that("dunnettTest: list and vector+g interfaces give identical results", {
  x <- c(2.9, 3.0, 2.5, 2.6, 3.2)
  y <- c(3.8, 2.7, 4.0, 2.4)
  z <- c(2.8, 3.4, 3.7, 2.2, 2.0)
  
  res_list <- dunnettTest(list(x, y, z))
  res_vec  <- dunnettTest(c(x, y, z), rep(1:3, c(5, 4, 5)))
  
  expect_equal(res_list[[1L]][, "diff"],  res_vec[[1L]][, "diff"],
               tolerance = 1e-10)
  expect_equal(res_list[[1L]][, "pval"],  res_vec[[1L]][, "pval"],
               tolerance = 1e-10)
})

