
# -------------------------------------------------------------------------
# Basic functionality
# -------------------------------------------------------------------------
test_that("pageTest returns an htest object", {
  
  soa.mat <- matrix(c(.797,.873,.888,.923,.942,.956,
                      .794,.772,.908,.982,.946,.913,
                      .838,.801,.853,.951,.883,.837,
                      .815,.801,.747,.859,.887,.902),
                    nrow = 4, byrow = TRUE)
  
  res <- pageTest(soa.mat)
  
  expect_s3_class(res, "htest")
  expect_named(res$statistic, "L")
  expect_named(res$parameter, c("k", "n"))
  expect_true(is.numeric(res$p.value))
  expect_gte(res$p.value, 0)
  expect_lte(res$p.value, 1)
})
# -------------------------------------------------------------------------
# Known reference values
# -------------------------------------------------------------------------
test_that("Siegel & Castellan example gives correct L and p-value", {
  
  # Craig's data, Siegel & Castellan p. 186
  soa.mat <- matrix(c(.797,.873,.888,.923,.942,.956,
                      .794,.772,.908,.982,.946,.913,
                      .838,.801,.853,.951,.883,.837,
                      .815,.801,.747,.859,.887,.902),
                    nrow = 4, byrow = TRUE)
  
  res <- pageTest(soa.mat)
  
  expect_equal(unname(res$statistic["L"]), 342)
  expect_equal(res$p.value, 0.000566074, tolerance = 1e-6)
  expect_match(res$method, "exact")
})

test_that("Sachs example gives correct L and p-value (k=4, n=5)", {
  
  # Sachs: k=4 groups, n=5 blocks, L=133, p > 0.05
  expect_equal(.pPage(k = 4L, n = 5L, L = 133L), 0.1265, tolerance = 1e-3)
})
# -------------------------------------------------------------------------
# Parameter output
# -------------------------------------------------------------------------
test_that("parameter k and n are correct", {
  
  soa.mat <- matrix(c(.797,.873,.888,.923,.942,.956,
                      .794,.772,.908,.982,.946,.913,
                      .838,.801,.853,.951,.883,.837,
                      .815,.801,.747,.859,.887,.902),
                    nrow = 4, byrow = TRUE)
  
  res <- pageTest(soa.mat)
  
  expect_equal(unname(res$parameter["k"]), 6L)
  expect_equal(unname(res$parameter["n"]), 4L)
})
# -------------------------------------------------------------------------
# Formula interface
# -------------------------------------------------------------------------
test_that("formula interface gives same result as default", {
  
  mat <- matrix(c(3,2,1,4, 4,2,3,1, 4,1,2,3, 4,2,3,1,
                  3,2,1,4, 4,1,2,3, 4,3,2,1, 3,1,2,4,
                  3,1,4,2),
                nrow = 9, byrow = TRUE,
                dimnames = list(1:9, LETTERS[1:4]))
  
  x <- mat[, c("B","C","D","A")]
  
  res_mat <- pageTest(x)
  
  plng <- data.frame(
    block = rep(1:9, 4),
    group = rep(c("B","C","D","A"), each = 9),
    y     = as.vector(x)
  )
  plng$group <- factor(plng$group, levels = c("B","C","D","A"))
  
  res_formula <- pageTest(y ~ group | block, data = plng)
  
  expect_equal(unname(res_mat$statistic), unname(res_formula$statistic))
  expect_equal(res_mat$p.value, res_formula$p.value)
})
# -------------------------------------------------------------------------
# Matrix vs vector interface
# -------------------------------------------------------------------------
test_that("matrix and vector interface give same result", {
  
  mat <- matrix(c(.797,.873,.888,.923,.942,.956,
                  .794,.772,.908,.982,.946,.913,
                  .838,.801,.853,.951,.883,.837,
                  .815,.801,.747,.859,.887,.902),
                nrow = 4, byrow = TRUE)
  
  res_mat <- pageTest(mat)
  
  res_vec <- pageTest(
    y      = as.vector(t(mat)),
    groups = factor(rep(seq_len(ncol(mat)), nrow(mat))),
    blocks = factor(rep(seq_len(nrow(mat)), each = ncol(mat)))
  )
  
  expect_equal(unname(res_mat$statistic), unname(res_vec$statistic))
  expect_equal(res_mat$p.value, res_vec$p.value)
})
# -------------------------------------------------------------------------
# NA handling
# -------------------------------------------------------------------------
test_that("incomplete blocks are removed", {
  
  mat <- matrix(c(.797,.873,.888,.923,.942,.956,
                  .794,.772,.908,.982,.946,.913,
                  .838,.801,.853,.951,.883,.837,
                  .815,.801,.747,.859,.887,.902),
                nrow = 4, byrow = TRUE)
  
  mat_na      <- mat
  mat_na[2,3] <- NA
  
  res_full <- pageTest(mat[-2, ])
  res_na   <- pageTest(mat_na)
  
  expect_equal(unname(res_full$statistic), unname(res_na$statistic))
  expect_equal(res_full$p.value, res_na$p.value)
})
# -------------------------------------------------------------------------
# Asymptotic method for k > 15
# -------------------------------------------------------------------------
test_that("asymptotic method used for k > 15", {
  
  set.seed(1)
  mat <- matrix(rnorm(16 * 5), nrow = 5, ncol = 16)
  
  res <- pageTest(mat)
  
  expect_match(res$method, "asymptotic")
  expect_gte(res$p.value, 0)
  expect_lte(res$p.value, 1)
})

test_that("exact method used for k <= 15", {
  
  soa.mat <- matrix(c(.797,.873,.888,.923,.942,.956,
                      .794,.772,.908,.982,.946,.913,
                      .838,.801,.853,.951,.883,.837,
                      .815,.801,.747,.859,.887,.902),
                    nrow = 4, byrow = TRUE)
  
  res <- pageTest(soa.mat)
  
  expect_match(res$method, "exact")
})
# -------------------------------------------------------------------------
# Detects increasing trend
# -------------------------------------------------------------------------
test_that("detects strong increasing trend", {
  
  set.seed(1)
  n <- 10; k <- 5
  mat <- matrix(0, nrow = n, ncol = k)
  for (i in seq_len(n))
    mat[i, ] <- rank(seq_len(k) + rnorm(k, sd = 0.1))
  
  res <- pageTest(mat)
  
  expect_lt(res$p.value, 0.01)
})

test_that("no trend gives non-significant p-value", {
  
  set.seed(42)
  mat <- matrix(rnorm(50), nrow = 10, ncol = 5)
  
  res <- pageTest(mat)
  
  expect_gt(res$p.value, 0.05)
})
# -------------------------------------------------------------------------
# Input validation
# -------------------------------------------------------------------------
test_that("k < 3 throws error", {
  
  mat <- matrix(rnorm(10), nrow = 5, ncol = 2)
  
  expect_error(pageTest(mat), "3")
})

test_that("length mismatch throws error", {
  
  expect_error(
    pageTest(y = 1:10, groups = 1:4, blocks = 1:10),
    "same length"
  )
})

test_that("replicated design throws error", {
  
  expect_error(
    pageTest(
      y      = 1:12,
      groups = rep(1:3, 4),
      blocks = rep(1:3, each = 4)
    ),
    "unreplicated"
  )
})
# -------------------------------------------------------------------------
# Print compatibility
# -------------------------------------------------------------------------
test_that("print.htest works", {
  
  soa.mat <- matrix(c(.797,.873,.888,.923,.942,.956,
                      .794,.772,.908,.982,.946,.913,
                      .838,.801,.853,.951,.883,.837,
                      .815,.801,.747,.859,.887,.902),
                    nrow = 4, byrow = TRUE)
  
  res <- pageTest(soa.mat)
  
  expect_output(print(res), "Page")
})

