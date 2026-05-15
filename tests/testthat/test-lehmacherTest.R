
# -------------------------------------------------------------------------
# Reference data
# -------------------------------------------------------------------------
x_ref <- matrix(c(400,40,20,10,
                  50,300,60,20,
                  10,40,120,5,
                  5,90,50,80), nrow = 4, byrow = TRUE)
# -------------------------------------------------------------------------
# Basic functionality
# -------------------------------------------------------------------------
test_that("lehmacherTest returns an mtest object", {
  
  res <- lehmacherTest(x_ref)
  
  expect_s3_class(res, "mtest")
  expect_true(is.numeric(res$statistic))
  expect_true(is.numeric(res$p.value))
  expect_true(is.numeric(res$p.value.corr))
  expect_equal(unname(res$parameter), 1L)
})
# -------------------------------------------------------------------------
# Reference values
# -------------------------------------------------------------------------
test_that("reference example gives correct chi-squared values", {
  
  res <- lehmacherTest(x_ref)
  
  expect_equal(unname(res$statistic),
               c(0.1852, 5.3333, 30.4054, 67.2222),
               tolerance = 1e-3)
})

test_that("reference example gives correct p-values", {
  
  res <- lehmacherTest(x_ref)
  
  expect_gt(res$p.value[1], 0.05)    # category 1: not significant
  expect_lt(res$p.value[2], 0.05)    # category 2: significant
  expect_lt(res$p.value[3], 0.001)   # category 3: highly significant
  expect_lt(res$p.value[4], 0.001)   # category 4: highly significant
})

test_that("adjusted p-values are >= raw p-values", {
  
  res <- lehmacherTest(x_ref)
  
  expect_true(all(res$p.value.corr >= res$p.value - 1e-10))
})
# -------------------------------------------------------------------------
# Perfect agreement gives statistic = 0
# -------------------------------------------------------------------------
test_that("perfect agreement gives statistic 0 for that category", {
  
  # diagonal = rowsum = colsum for category 1
  mat <- matrix(c(10,0,0,0,
                  0,5,2,1,
                  0,1,8,0,
                  0,2,1,6), nrow = 4, byrow = TRUE)
  
  res <- lehmacherTest(mat)
  
  expect_equal(unname(res$statistic[1]), 0)
})
# -------------------------------------------------------------------------
# Vector interface
# -------------------------------------------------------------------------
test_that("matrix and vector interface give same result", {
  
  d <- data.frame(
    expand.grid(r = c("A","B","C"), c = c("A","B","C"))[
      rep(1:9, times = c(20,3,0,10,30,5,5,15,40)), ]
  )
  
  res_mat <- lehmacherTest(
    matrix(c(20,3,0,10,30,5,5,15,40), nrow=3, byrow=TRUE)
  )
  res_vec <- lehmacherTest(d$r, d$c)
  
  expect_equal(unname(res_mat$statistic), unname(res_vec$statistic),
               tolerance = 1e-6)
})
# -------------------------------------------------------------------------
# Level unification
# -------------------------------------------------------------------------
test_that("mismatched but overlapping levels are unified", {
  
  x <- factor(c("A","A","B","B"), levels = c("A","B"))
  y <- factor(c("A","B","B","C"), levels = c("B","C"))
  
  expect_no_error(lehmacherTest(x, y))
})
# -------------------------------------------------------------------------
# Input validation
# -------------------------------------------------------------------------
test_that("non-square matrix throws error", {
  expect_error(lehmacherTest(matrix(1:6, nrow=2)), "square")
})

test_that("negative entries throw error", {
  mat      <- x_ref
  mat[1,1] <- -1
  expect_error(lehmacherTest(mat), "nonnegative")
})

test_that("non-integer counts produce warning", {
  mat      <- x_ref
  mat[1,1] <- 400.5
  expect_warning(lehmacherTest(mat), "non-integer")
})

test_that("missing y throws error", {
  expect_error(lehmacherTest(factor(c("A","B","A"))), "'y' must be given")
})

test_that("different length x and y throws error", {
  expect_error(
    lehmacherTest(factor(c("A","B")), factor(c("A","B","C"))),
    "same length"
  )
})
# -------------------------------------------------------------------------
# Print method
# -------------------------------------------------------------------------
test_that("print.mtest works", {
  
  res <- lehmacherTest(x_ref)
  
  expect_output(print(res), "Lehmacher")
  expect_output(print(res), "X-squared")
  expect_output(print(res), "pval")
})