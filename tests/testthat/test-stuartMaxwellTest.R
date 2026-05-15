

# -------------------------------------------------------------------------
# Basic functionality
# -------------------------------------------------------------------------
test_that("stuartMaxwellTest returns an htest object", {
  
  hyp <- as.table(matrix(c(20,3,0,10,30,5,5,15,40), nrow = 3))
  
  res <- stuartMaxwellTest(hyp)
  
  expect_s3_class(res, "htest")
  expect_named(res$statistic, "chi-squared")
  expect_named(res$parameter, "df")
  expect_true(is.numeric(res$p.value))
  expect_gte(res$p.value, 0)
  expect_lte(res$p.value, 1)
})
# -------------------------------------------------------------------------
# Reference values
# -------------------------------------------------------------------------
test_that("Uebersax example gives correct chi-squared and p-value", {
  
  hyp <- as.table(matrix(c(20,3,0,10,30,5,5,15,40), nrow = 3))
  
  res <- stuartMaxwellTest(hyp)
  
  expect_equal(unname(res$statistic["chi-squared"]), 13.765, tolerance = 1e-3)
  expect_equal(res$p.value, 0.001025728, tolerance = 1e-5)
  expect_equal(unname(res$parameter["df"]), 2L)
})

test_that("Agresti 4x4 example gives correct chi-squared and p-value", {
  
  mc <- as.table(matrix(c(
    732,1524,1575,1577,1602,837,1554,1437,
    1672,1600,841,1363,1385,1484,1524,791), nrow = 4))
  
  res <- stuartMaxwellTest(mc)
  
  expect_equal(unname(res$statistic["chi-squared"]), 0.089722, tolerance = 1e-4)
  expect_equal(res$p.value, 0.993, tolerance = 1e-2)
  expect_equal(unname(res$parameter["df"]), 3L)
})
# -------------------------------------------------------------------------
# Matrix and vector interface give same result
# -------------------------------------------------------------------------
test_that("matrix and vector interface give same result", {
  
  hyp <- as.table(matrix(c(20,3,0,10,30,5,5,15,40), nrow = 3))
  
  d.hyp <- data.frame(
    expand.grid(r = c("A","B","C"), c = c("A","B","C"))[
      rep(1:9, times = as.vector(hyp)), ]
  )
  
  res_mat <- stuartMaxwellTest(hyp)
  res_vec <- stuartMaxwellTest(d.hyp$r, d.hyp$c)
  
  expect_equal(unname(res_mat$statistic), unname(res_vec$statistic),
               tolerance = 1e-6)
  expect_equal(res_mat$p.value, res_vec$p.value, tolerance = 1e-6)
})
# -------------------------------------------------------------------------
# Level unification
# -------------------------------------------------------------------------
test_that("mismatched but overlapping levels are unified", {
  
  # x has A,B; y has B,C -> union = A,B,C
  x <- factor(c("A","A","B","B"), levels = c("A","B"))
  y <- factor(c("A","B","B","C"), levels = c("B","C"))
  
  # Should not throw error about level mismatch
  expect_no_error(stuartMaxwellTest(x, y))
})

test_that("same levels give same result as explicit matrix", {
  
  x <- c("A","A","B","C","A","B")
  y <- c("A","B","A","C","B","A")
  
  res_vec <- stuartMaxwellTest(factor(x), factor(y))
  res_mat <- stuartMaxwellTest(table(factor(x), factor(y)))
  
  expect_equal(unname(res_vec$statistic), unname(res_mat$statistic),
               tolerance = 1e-6)
})
# -------------------------------------------------------------------------
# Perfect agreement handling
# -------------------------------------------------------------------------
test_that("perfect agreement in one category reduces matrix but df stays k-1", {
  
  # Category C has perfect agreement (diagonal = row = col total)
  mat <- matrix(c(10,3,0, 2,8,0, 0,0,15), nrow = 3,
                dimnames = list(c("A","B","C"), c("A","B","C")))
  
  res <- stuartMaxwellTest(mat)
  
  # df should still be k_original - 1 = 2
  expect_equal(unname(res$parameter["df"]), 2L)
  expect_true(is.finite(res$p.value))
})

test_that("too many perfect agreement categories throws error", {
  
  # All three categories have perfect agreement
  mat <- diag(c(10, 20, 30))
  
  expect_error(stuartMaxwellTest(mat), "perfect agreement")
})
# -------------------------------------------------------------------------
# Input validation
# -------------------------------------------------------------------------
test_that("non-square matrix throws error", {
  
  mat <- matrix(1:6, nrow = 2)
  
  expect_error(stuartMaxwellTest(mat), "square")
})

test_that("non-numeric matrix throws error", {
  
  mat <- matrix(c("1","2","3","4"), nrow = 2)
  
  expect_error(stuartMaxwellTest(mat), "numeric")
})

test_that("negative entries throw error", {
  
  mat <- matrix(c(10,-1,2,5), nrow = 2)
  
  expect_error(stuartMaxwellTest(mat), "nonnegative")
})

test_that("missing y throws error", {
  
  expect_error(stuartMaxwellTest(factor(c("A","B","A"))), "'y' must be given")
})

test_that("different length x and y throws error", {
  
  expect_error(
    stuartMaxwellTest(factor(c("A","B")), factor(c("A","B","C"))),
    "same length"
  )
})

test_that("fewer than 2 levels throws error", {
  
  expect_error(
    stuartMaxwellTest(factor(c("A","A")), factor(c("A","A"))),
    "2 distinct"
  )
})

test_that("singular matrix throws error", {
  
  # Perfectly correlated -> singular S
  mat <- matrix(c(10,0,0, 0,10,0, 0,0,10), nrow = 3)
  
  expect_error(stuartMaxwellTest(mat), "singular|perfect agreement")
})

test_that("mismatched dimnames produce warning", {
  
  mat <- matrix(c(10,3,2,8), nrow = 2,
                dimnames = list(c("A","B"), c("X","Y")))
  
  expect_warning(stuartMaxwellTest(mat), "differ")
})
# -------------------------------------------------------------------------
# n in output
# -------------------------------------------------------------------------
test_that("n equals total sample size", {
  
  hyp <- as.table(matrix(c(20,3,0,10,30,5,5,15,40), nrow = 3))
  
  res <- stuartMaxwellTest(hyp)
  
  expect_equal(res$n, sum(hyp))
})
# -------------------------------------------------------------------------
# 2x2 equals McNemar
# -------------------------------------------------------------------------
test_that("2x2 case gives same p-value as mcnemar.test", {
  
  mat <- matrix(c(10,3,7,20), nrow = 2)
  
  res_sm  <- stuartMaxwellTest(mat)
  res_mc  <- mcnemar.test(mat, correct = FALSE)
  
  expect_equal(res_sm$p.value, res_mc$p.value, tolerance = 1e-6)
})
# -------------------------------------------------------------------------
# Print compatibility
# -------------------------------------------------------------------------
test_that("print.htest works", {
  
  hyp <- as.table(matrix(c(20,3,0,10,30,5,5,15,40), nrow = 3))
  
  res <- stuartMaxwellTest(hyp)
  
  expect_output(print(res), "Stuart-Maxwell")
})

