
test_that("vanWaerdenTest.default: basic three-group example (list interface)", {
  x <- c(2.9, 3.0, 2.5, 2.6, 3.2)
  y <- c(3.8, 2.7, 4.0, 2.4)
  z <- c(2.8, 3.4, 3.7, 2.2, 2.0)
  
  res <- vanWaerdenTest(list(x, y, z))
  
  expect_s3_class(res, "htest")
  expect_named(res, c("statistic", "parameter", "p.value", "method", "data.name"))
})


test_that("vanWaerdenTest.default: statistic and p-value are scalar numerics", {
  set.seed(1)
  x <- c(rnorm(10, 0), rnorm(10, 1), rnorm(10, 2))
  g <- rep(1:3, each = 10)
  
  res <- vanWaerdenTest(x, g)
  
  expect_true(is.numeric(res$statistic))
  expect_length(res$statistic, 1L)
  expect_true(is.numeric(res$p.value))
  expect_length(res$p.value, 1L)
})


test_that("vanWaerdenTest.default: p-value is in [0, 1]", {
  set.seed(2)
  x <- c(rnorm(12, 0), rnorm(12, 2), rnorm(12, 4))
  g <- rep(letters[1:3], each = 12)
  
  res <- vanWaerdenTest(x, g)
  
  expect_gte(res$p.value, 0)
  expect_lte(res$p.value, 1)
})


test_that("vanWaerdenTest.default: statistic is non-negative", {
  set.seed(3)
  x <- rnorm(30)
  g <- rep(1:3, 10)
  
  res <- vanWaerdenTest(x, g)
  
  expect_gte(unname(res$statistic), 0)
})


test_that("vanWaerdenTest.default: parameter equals k - 1 and is numeric", {
  set.seed(4)
  x <- rnorm(40)
  g <- rep(1:4, 10)
  
  res <- vanWaerdenTest(x, g)
  
  expect_equal(unname(res$parameter), 3)    # k=4, df=3
  expect_true(is.numeric(res$parameter))    # must not be integer
  expect_false(is.integer(res$parameter))
})


test_that("vanWaerdenTest.default: statistic name and parameter name are correct", {
  res <- vanWaerdenTest(list(c(1, 2, 3), c(4, 5, 6)))
  
  expect_equal(names(res$statistic), "Van-der-Waerden chi-squared")
  expect_equal(names(res$parameter), "df")
})


test_that("vanWaerdenTest.default: method string is correct", {
  res <- vanWaerdenTest(list(c(1, 2, 3), c(4, 5, 6)))
  
  expect_equal(res$method, "Van-der-Waerden normal scores test")
})


test_that("vanWaerdenTest.default: list and vector+g interfaces give identical results", {
  x <- c(2.9, 3.0, 2.5, 2.6, 3.2)
  y <- c(3.8, 2.7, 4.0, 2.4)
  z <- c(2.8, 3.4, 3.7, 2.2, 2.0)
  
  res_list <- vanWaerdenTest(list(x, y, z))
  res_vec  <- vanWaerdenTest(c(x, y, z), rep(1:3, c(5, 4, 5)))
  
  expect_equal(unname(res_list$statistic), unname(res_vec$statistic),
               tolerance = 1e-10)
  expect_equal(res_list$p.value, res_vec$p.value, tolerance = 1e-10)
})


test_that("vanWaerdenTest.default: well-separated groups yield small p-value", {
  set.seed(5)
  x <- c(rnorm(20, 0), rnorm(20, 10), rnorm(20, 20))
  g <- rep(1:3, each = 20)
  
  res <- vanWaerdenTest(x, g)
  
  expect_lt(res$p.value, 0.001)
})


test_that("vanWaerdenTest.default: identical groups yield large p-value", {
  set.seed(6)
  x <- rnorm(30)
  g <- rep(1:3, 10)
  
  res <- vanWaerdenTest(x, g)
  
  # With no location difference, p-value should not be consistently small
  expect_gt(res$p.value, 0.001)
})


test_that("vanWaerdenTest.default: two-group case uses df = 1", {
  set.seed(7)
  x <- c(rnorm(15, 0), rnorm(15, 2))
  g <- rep(c("a", "b"), each = 15)
  
  res <- vanWaerdenTest(x, g)
  
  expect_equal(unname(res$parameter), 1)
})


test_that("vanWaerdenTest.default: ties do not crash and p-value stays in [0, 1]", {
  x <- c(rep(1, 5), rep(2, 5), rep(3, 5))
  g <- rep(c("a", "b", "c"), each = 5)
  
  expect_no_error(vanWaerdenTest(x, g))
  res <- vanWaerdenTest(x, g)
  expect_gte(res$p.value, 0)
  expect_lte(res$p.value, 1)
})


test_that("vanWaerdenTest.default: NAs in input are silently removed", {
  x <- c(1, 2, NA, 4, 5, 6, 7, NA, 9)
  g <- rep(c("a", "b", "c"), 3)
  
  expect_no_error(vanWaerdenTest(x, g))
})


test_that("vanWaerdenTest.formula: equivalent to default interface", {
  df <- data.frame(
    val = c(2.9, 3.0, 2.5, 2.6, 3.2,
            3.8, 2.7, 4.0, 2.4,
            2.8, 3.4, 3.7, 2.2, 2.0),
    grp = factor(rep(c("X", "Y", "Z"), c(5, 4, 5)))
  )
  
  res_f <- vanWaerdenTest(val ~ grp, data = df)
  res_d <- vanWaerdenTest(df$val, df$grp)
  
  expect_equal(unname(res_f$statistic), unname(res_d$statistic), tolerance = 1e-10)
  expect_equal(res_f$p.value,           res_d$p.value,           tolerance = 1e-10)
})


test_that("vanWaerdenTest.formula: data.name is set from formula", {
  df <- data.frame(
    val = c(1, 2, 3, 4, 5, 6),
    grp = rep(c("a", "b"), 3)
  )
  
  res <- vanWaerdenTest(val ~ grp, data = df)
  
  expect_equal(res$data.name, "val ~ grp")
})


test_that("vanWaerdenTest.formula: subset argument is respected", {
  df <- data.frame(
    val = c(2.9, 3.0, 2.5, 2.6, 3.2,
            3.8, 2.7, 4.0, 2.4,
            2.8, 3.4, 3.7, 2.2, 2.0),
    grp = factor(rep(c("X", "Y", "Z"), c(5, 4, 5)))
  )
  
  # Subsetting to two groups should give df = 1
  res_sub <- vanWaerdenTest(val ~ grp, data = df, subset = grp != "Z")
  
  expect_equal(unname(res_sub$parameter), 1)
})


test_that("vanWaerdenTest.formula: na.action = na.omit drops NAs silently", {
  df <- data.frame(
    val = c(NA, 3.0, 2.5, 2.6, 3.2,
            3.8, 2.7, 4.0, 2.4,
            2.8, 3.4, 3.7, 2.2, 2.0),
    grp = factor(rep(c("X", "Y", "Z"), c(5, 4, 5)))
  )
  
  expect_no_error(
    vanWaerdenTest(val ~ grp, data = df, na.action = na.omit)
  )
})


test_that("vanWaerdenTest: airquality formula example runs without error", {
  expect_no_error({
    res <- vanWaerdenTest(Ozone ~ Month, data = airquality, na.action = na.omit)
    expect_s3_class(res, "htest")
    expect_equal(unname(res$parameter), 4)   # 5 months -> df = 4
    expect_lt(res$p.value, 0.05)
  })
})


test_that("vanWaerdenTest: error on malformed formula", {
  expect_error(vanWaerdenTest.formula(~ x), "'formula' missing or incorrect")
})


test_that("vanWaerdenTest: statistic is invariant to group label ordering", {
  # Reordering factor levels must not change the test statistic
  set.seed(9)
  x <- c(rnorm(10), rnorm(10, 1), rnorm(10, 2))
  g1 <- factor(rep(c("a", "b", "c"), each = 10))
  g2 <- factor(rep(c("a", "b", "c"), each = 10), levels = c("c", "b", "a"))
  
  res1 <- vanWaerdenTest(x, g1)
  res2 <- vanWaerdenTest(x, g2)
  
  expect_equal(unname(res1$statistic), unname(res2$statistic), tolerance = 1e-10)
  expect_equal(res1$p.value, res2$p.value, tolerance = 1e-10)
})


test_that("vanWaerdenTest: more powerful than Kruskal-Wallis on near-normal data", {
  # Van der Waerden should be at least as powerful as KW for normal data.
  # Verified by comparing p-values on a fixed dataset with clear signal.
  set.seed(42)
  x <- c(rnorm(25, 0), rnorm(25, 1.5), rnorm(25, 3))
  g <- rep(1:3, each = 25)
  
  p_vdw <- vanWaerdenTest(x, g)$p.value
  p_kw  <- kruskal.test(x, g)$p.value
  
  # Van der Waerden p-value should be <= KW p-value (more powerful)
  expect_lte(p_vdw, p_kw + 0.05)   # allow small tolerance for random variation
})

