

# -------------------------------------------------------------------------
# Basic functionality
# -------------------------------------------------------------------------
test_that("bhapkarTest returns an htest object", {
  
  mc <- as.table(matrix(c(20,3,0,10,30,5,5,15,40), nrow = 3))
  
  res <- bhapkarTest(mc)
  
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
  
  mc <- as.table(matrix(c(20,3,0,10,30,5,5,15,40), nrow = 3))
  
  res <- bhapkarTest(mc)
  
  expect_equal(unname(res$statistic["chi-squared"]), 15.423, tolerance = 1e-3)
  expect_equal(res$p.value, 0.000447588, tolerance = 1e-6)
  expect_equal(unname(res$parameter["df"]), 2L)
})
# -------------------------------------------------------------------------
# Bhapkar >= Stuart-Maxwell
# -------------------------------------------------------------------------
test_that("Bhapkar statistic is larger than Stuart-Maxwell", {
  
  mc <- as.table(matrix(c(20,3,0,10,30,5,5,15,40), nrow = 3))
  
  res_sm <- stuartMaxwellTest(mc)
  res_bh <- bhapkarTest(mc)
  
  expect_gt(
    unname(res_bh$statistic["chi-squared"]),
    unname(res_sm$statistic["chi-squared"])
  )
})

test_that("Bhapkar p-value <= Stuart-Maxwell p-value", {
  
  mc <- as.table(matrix(c(20,3,0,10,30,5,5,15,40), nrow = 3))
  
  res_sm <- stuartMaxwellTest(mc)
  res_bh <- bhapkarTest(mc)
  
  expect_lte(res_bh$p.value, res_sm$p.value)
})
# -------------------------------------------------------------------------
# Asymptotic equivalence for large n
# -------------------------------------------------------------------------
test_that("Bhapkar and Stuart-Maxwell p-values converge for large n", {
  
  mc_large <- as.table(matrix(c(20,3,0,10,30,5,5,15,40) * 1000L, nrow = 3))
  
  res_sm <- stuartMaxwellTest(mc_large)
  res_bh <- bhapkarTest(mc_large)
  
  # Both should be essentially 0 for large n with clear signal
  expect_lt(res_bh$p.value, 1e-10)
  expect_lt(res_sm$p.value, 1e-10)
})

# -------------------------------------------------------------------------
# Matrix and vector interface
# -------------------------------------------------------------------------
test_that("matrix and vector interface give same result", {
  
  mc  <- as.table(matrix(c(20,3,0,10,30,5,5,15,40), nrow = 3))
  
  d <- data.frame(
    expand.grid(r = c("A","B","C"), c = c("A","B","C"))[
      rep(1:9, times = as.vector(mc)), ]
  )
  
  res_mat <- bhapkarTest(mc)
  res_vec <- bhapkarTest(d$r, d$c)
  
  expect_equal(unname(res_mat$statistic), unname(res_vec$statistic),
               tolerance = 1e-6)
  expect_equal(res_mat$p.value, res_vec$p.value, tolerance = 1e-6)
})
# -------------------------------------------------------------------------
# df equals k - 1
# -------------------------------------------------------------------------
test_that("df equals k - 1", {
  
  mc <- as.table(matrix(c(20,3,0,10,30,5,5,15,40), nrow = 3))
  
  res <- bhapkarTest(mc)
  
  expect_equal(unname(res$parameter["df"]), 2L)
})
# -------------------------------------------------------------------------
# 2x2 case
# -------------------------------------------------------------------------
test_that("2x2 gives finite p-value", {
  
  mat <- matrix(c(10,3,7,20), nrow = 2)
  
  res <- bhapkarTest(mat)
  
  expect_true(is.finite(res$p.value))
  expect_lt(res$p.value, 1)
})
# -------------------------------------------------------------------------
# method string
# -------------------------------------------------------------------------
test_that("method string contains Bhapkar", {
  
  mc <- as.table(matrix(c(20,3,0,10,30,5,5,15,40), nrow = 3))
  
  res <- bhapkarTest(mc)
  
  expect_match(res$method, "Bhapkar")
})
# -------------------------------------------------------------------------
# Print compatibility
# -------------------------------------------------------------------------
test_that("print.htest works", {
  
  mc <- as.table(matrix(c(20,3,0,10,30,5,5,15,40), nrow = 3))
  
  res <- bhapkarTest(mc)
  
  expect_output(print(res), "Bhapkar")
})