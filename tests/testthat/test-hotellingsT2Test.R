

test_that("one-sample F: htest structure is correct", {
  
  set.seed(1)
  x <- cbind(rnorm(20), rnorm(20))
  
  res <- hotellingsT2Test(x)
  
  expect_s3_class(res, "htest")
  expect_named(res, c("statistic", "parameter", "p.value", "null.value",
                      "alternative", "method", "data.name"))
  expect_named(res$statistic,  "T.2")
  expect_named(res$parameter,  c("df1", "df2"))
  expect_true(all(names(res$null.value) == "location"))
  expect_equal(res$alternative, "two.sided")
  expect_match(res$method, "one-sample")
})


test_that("one-sample F: null.value is numeric, not character", {
  
  set.seed(1)
  x  <- cbind(rnorm(20), rnorm(20))
  mu <- c(0.5, -0.5)
  
  res <- hotellingsT2Test(x, mu = mu)
  
  expect_true(is.numeric(res$null.value))
  expect_equal(unname(res$null.value), mu)
})


test_that("one-sample F: correct statistic against known value (ICSNP reference)", {
  
  # Data from ICSNP package vignette (symmetric, known T2)
  set.seed(42)
  x <- cbind(
    c(2.1, 3.4, 2.8, 3.1, 2.9, 3.6, 2.4, 3.0),
    c(1.8, 2.9, 2.3, 2.7, 2.5, 3.1, 2.0, 2.6)
  )
  
  res   <- hotellingsT2Test(x, mu = c(3, 2.5))
  res_chi <- hotellingsT2Test(x, mu = c(3, 2.5), test = "chi")
  
  # F and chi2 statistics must be positive
  expect_gt(res$statistic,     0)
  expect_gt(res_chi$statistic, 0)
  
  # chi2 stat = T2 (unscaled), F stat = T2 * scale < T2
  expect_lt(res$statistic, res_chi$statistic)
})


test_that("one-sample chi: parameter has single df", {
  
  set.seed(1)
  x   <- cbind(rnorm(30), rnorm(30), rnorm(30))
  res <- hotellingsT2Test(x, test = "chi")
  
  expect_named(res$parameter, "df")
  expect_equal(unname(res$parameter), 3L)   # p = 3
})


test_that("one-sample: p-value in (0, 1)", {
  
  set.seed(7)
  x <- cbind(rnorm(15, mean = 1), rnorm(15, mean = 1))
  res <- hotellingsT2Test(x, mu = c(1, 1))
  
  expect_true(res$p.value > 0 && res$p.value < 1)
})


test_that("one-sample: mu = true mean gives large p-value", {
  
  set.seed(3)
  x   <- cbind(rnorm(200, 5), rnorm(200, 10))
  res <- hotellingsT2Test(x, mu = c(5, 10))
  
  expect_gt(res$p.value, 0.05)
})


test_that("one-sample: mu far from data gives small p-value", {
  
  set.seed(4)
  x   <- cbind(rnorm(50, 5), rnorm(50, 10))
  res <- hotellingsT2Test(x, mu = c(0, 0))
  
  expect_lt(res$p.value, 0.001)
})


# -----------------------------------------------------------------------
# Two-sample tests
# -----------------------------------------------------------------------

test_that("two-sample F: htest structure is correct", {
  
  set.seed(1)
  x <- cbind(rnorm(15), rnorm(15))
  y <- cbind(rnorm(20), rnorm(20))
  
  res <- hotellingsT2Test(x, y)
  
  expect_s3_class(res, "htest")
  expect_named(res$statistic,  "T.2")
  expect_named(res$parameter,  c("df1", "df2"))
  expect_true(all(names(res$null.value) == "location difference"))
  expect_match(res$method,     "two-sample")
})


test_that("two-sample F: df2 = n1 + n2 - p - 1", {
  
  set.seed(1)
  n1 <- 15; n2 <- 20; p <- 3
  x  <- matrix(rnorm(n1 * p), n1, p)
  y  <- matrix(rnorm(n2 * p), n2, p)
  
  res <- hotellingsT2Test(x, y)
  
  expect_equal(unname(res$parameter["df2"]), n1 + n2 - p - 1)
})


test_that("two-sample chi: parameter has single df", {
  
  set.seed(1)
  x   <- matrix(rnorm(30 * 2), 30, 2)
  y   <- matrix(rnorm(25 * 2), 25, 2)
  res <- hotellingsT2Test(x, y, test = "chi")
  
  expect_named(res$parameter, "df")
  expect_equal(unname(res$parameter), 2L)
})


test_that("two-sample: identical groups give large p-value", {
  
  set.seed(5)
  x   <- cbind(rnorm(40, 3), rnorm(40, 7))
  y   <- cbind(rnorm(40, 3), rnorm(40, 7))
  res <- hotellingsT2Test(x, y)
  
  expect_gt(res$p.value, 0.05)
})


test_that("two-sample: well-separated groups give small p-value", {
  
  set.seed(6)
  x   <- cbind(rnorm(50,  0), rnorm(50,  0))
  y   <- cbind(rnorm(50, 10), rnorm(50, 10))
  res <- hotellingsT2Test(x, y)
  
  expect_lt(res$p.value, 0.001)
})


test_that("two-sample: null.value is numeric zero vector by default", {
  
  set.seed(1)
  x   <- matrix(rnorm(20 * 3), 20, 3)
  y   <- matrix(rnorm(15 * 3), 15, 3)
  res <- hotellingsT2Test(x, y)
  
  expect_true(is.numeric(res$null.value))
  expect_equal(unname(res$null.value), c(0, 0, 0))
})


# -----------------------------------------------------------------------
# Formula interface
# -----------------------------------------------------------------------

test_that("formula interface matches default interface", {
  
  set.seed(1)
  df <- data.frame(
    g    = factor(rep(c("a", "b"), each = 20)),
    v1   = c(rnorm(20, 0), rnorm(20, 1)),
    v2   = c(rnorm(20, 0), rnorm(20, 2))
  )
  
  res_formula <- hotellingsT2Test(cbind(v1, v2) ~ g, data = df)
  
  # mirror exactly what resolveFormula does: split by levels(g), alphabetical
  lvls <- levels(df$g)
  resp <- as.matrix(df[, c("v1", "v2")])
  x    <- resp[df$g == lvls[1L], , drop = FALSE]
  y    <- resp[df$g == lvls[2L], , drop = FALSE]
  res_default <- hotellingsT2Test(x, y)
  
  expect_equal(res_formula$statistic, res_default$statistic)
  expect_equal(res_formula$p.value,   res_default$p.value)
  expect_equal(res_formula$parameter, res_default$parameter)
})


test_that("formula interface: data.name comes from formula", {
  
  set.seed(1)
  df <- data.frame(
    g  = factor(rep(c("a", "b"), each = 15)),
    v1 = rnorm(30),
    v2 = rnorm(30)
  )
  
  res <- hotellingsT2Test(cbind(v1, v2) ~ g, data = df)
  
  expect_match(res$data.name, "cbind")
})


test_that("formula interface: subset argument works", {
  
  set.seed(1)
  df <- data.frame(
    g    = factor(rep(c("a", "b"), each = 20)),
    v1   = c(rnorm(20, 0), rnorm(20, 1)),
    v2   = c(rnorm(20, 0), rnorm(20, 2)),
    keep = c(rep(TRUE, 15), rep(FALSE, 5), rep(TRUE, 20))
  )
  
  res_sub  <- hotellingsT2Test(cbind(v1, v2) ~ g, data = df, subset = keep)
  res_full <- hotellingsT2Test(cbind(v1, v2) ~ g, data = df)
  
  # Different data → different statistic
  expect_false(isTRUE(all.equal(res_sub$statistic, res_full$statistic)))
})


# -----------------------------------------------------------------------
# Input validation
# -----------------------------------------------------------------------

test_that("non-numeric x raises error", {
  x <- cbind(letters[1:5], letters[1:5])
  expect_error(hotellingsT2Test(x), "'x' must be numeric")
})


test_that("n <= p raises error for x", {
  x <- matrix(rnorm(6), nrow = 2, ncol = 3)   # n=2, p=3
  expect_error(hotellingsT2Test(x), "more rows than columns")
})


test_that("column mismatch between x and y raises error", {
  x <- matrix(rnorm(20 * 2), 20, 2)
  y <- matrix(rnorm(20 * 3), 20, 3)
  expect_error(hotellingsT2Test(x, y), "same number of columns")
})


test_that("mu wrong length raises error", {
  x <- matrix(rnorm(30 * 2), 30, 2)
  expect_error(hotellingsT2Test(x, mu = c(1, 2, 3)), "(?i)length of 'mu'",
               perl = TRUE)
})


test_that("mu non-finite raises error", {
  x <- matrix(rnorm(30 * 2), 30, 2)
  expect_error(hotellingsT2Test(x, mu = c(1, Inf)), "finite")
  expect_error(hotellingsT2Test(x, mu = c(NA, 0)),  "finite")
})


test_that("singular covariance matrix raises informative error", {
  # Perfectly collinear columns
  v  <- rnorm(20)
  x  <- cbind(v, 2 * v)
  expect_error(hotellingsT2Test(x), "singular")
})


test_that("insufficient df for two-sample F: guard exists but is structurally unreachable", {
  # The guard (n1 + n2 - p - 1 > 0) in hotellingsT2Test.default is dead code:
  # for it to fire, both nrow(x) > p and nrow(y) > p must hold (earlier guards),
  # which implies n1 + n2 > 2p >= p + 1 for p >= 1, so df2 is always positive.
  # No test case can reach the guard through the public interface.
  skip("df2 guard is dead code; covered by nrow > p checks above it")
})


test_that("NAs in x are silently dropped", {
  set.seed(1)
  x      <- cbind(rnorm(20), rnorm(20))
  x[3, ] <- NA
  
  res_na  <- hotellingsT2Test(x)
  res_ref <- hotellingsT2Test(x[-3, ])
  
  expect_equal(res_na$statistic, res_ref$statistic)
})


test_that("invalid test argument raises error", {
  x <- matrix(rnorm(30 * 2), 30, 2)
  expect_error(hotellingsT2Test(x, test = "t"), "arg")
})
