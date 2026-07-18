
Job <- matrix(c(1,2,1,0, 3,3,6,1, 10,10,14,9, 6,7,12,11), 4, 4,
              dimnames = list(
                income       = c("< 15k","15-25k","25-40k","> 40k"),
                satisfaction = c("VeryD","LittleD","ModerateS","VeryS")))
# -------------------------------------------------------------------------
# Basic functionality
# -------------------------------------------------------------------------
test_that("mantelTrendTest returns an htest object", {
  
  res <- mantelTrendTest(Job)
  
  expect_s3_class(res, "htest")
  expect_named(res$statistic, "X-squared")
  expect_named(res$parameter, "df")
  expect_equal(unname(res$parameter), 1L)
  expect_true(is.numeric(res$p.value))
  expect_gte(res$p.value, 0)
  expect_lte(res$p.value, 1)
})
# -------------------------------------------------------------------------
# Reference values
# -------------------------------------------------------------------------
test_that("default scores give correct result", {
  
  res <- mantelTrendTest(Job)
  
  expect_equal(unname(res$statistic), 2.983, tolerance = 1e-3)
  expect_equal(res$p.value, 0.08414, tolerance = 1e-4)
})

test_that("custom srow scores give correct result", {
  
  res <- mantelTrendTest(Job, srow = c(7.5,20,32.5,60))
  
  expect_equal(unname(res$statistic), 3.8075, tolerance = 1e-3)
  expect_equal(res$p.value, 0.05102, tolerance = 1e-4)
})

test_that("custom scores change result vs default", {
  
  res_default <- mantelTrendTest(Job)
  res_custom  <- mantelTrendTest(Job, srow = c(7.5,20,32.5,60))
  
  expect_false(
    isTRUE(all.equal(
      unname(res_default$statistic),
      unname(res_custom$statistic)
    ))
  )
})
# -------------------------------------------------------------------------
# Pearson correlation
# -------------------------------------------------------------------------
test_that(".pearsonCor returns value in [-1, 1]", {
  
  r <- .pearsonCor(Job)
  
  expect_gte(r, -1)
  expect_lte(r, 1)
})

test_that(".pearsonCor changes with custom scores", {
  
  r_default <- .pearsonCor(Job)
  r_custom  <- .pearsonCor(Job, srow = c(7.5,20,32.5,60))
  
  expect_false(isTRUE(all.equal(r_default, r_custom)))
})

test_that("Q = (n-1) * r^2", {
  
  r   <- .pearsonCor(Job)
  res <- mantelTrendTest(Job)
  
  expect_equal(
    unname(res$statistic),
    (sum(Job) - 1) * r^2,
    tolerance = 1e-10
  )
})
# -------------------------------------------------------------------------
# Score length validation
# -------------------------------------------------------------------------
test_that("wrong srow length throws error", {
  expect_error(mantelTrendTest(Job, srow = 1:3), "srow")
})

test_that("wrong scol length throws error", {
  expect_error(mantelTrendTest(Job, scol = 1:3), "scol")
})
# -------------------------------------------------------------------------
# Input validation
# -------------------------------------------------------------------------
test_that("non-numeric matrix throws error", {
  expect_error(mantelTrendTest(matrix(c("a","b","c","d"), 2)), "numeric")
})

test_that("negative entries throw error", {
  mat      <- Job
  mat[1,1] <- -1
  expect_error(mantelTrendTest(mat), "nonnegative")
})
# -------------------------------------------------------------------------
# Print compatibility
# -------------------------------------------------------------------------
test_that("print.htest works", {
  res <- mantelTrendTest(Job)
  expect_output(print(res), "Mantel")
})


test_that("method string reflects linear association", {
  res <- mantelTrendTest(Job)
  expect_match(res$method, "linear")
})

test_that("estimate contains Pearson r", {
  res <- mantelTrendTest(Job)
  expect_named(res$estimate, "r")
  expect_gte(unname(res$estimate), -1)
  expect_lte(unname(res$estimate),  1)
})

test_that("Q = (n-1) * r^2", {
  res <- mantelTrendTest(Job)
  expect_equal(
    unname(res$statistic),
    (sum(Job) - 1) * unname(res$estimate)^2,
    tolerance = 1e-10
  )
})

test_that("empty table throws error", {
  expect_error(mantelTrendTest(matrix(0, 3, 3)), "2 observations")
})

test_that("zero variance scores throw error", {
  expect_error(
    mantelTrendTest(Job, srow = rep(1, 4)),
    "zero variance"
  )
})

test_that("non-monotone scores produce warning", {
  # c(4,3,2,1) is strictly decreasing, which is monotone (r^2 is
  # sign-invariant) and must NOT warn; a genuinely non-monotone sequence
  # is needed to exercise the warning path
  expect_warning(
    mantelTrendTest(Job, srow = c(1,3,2,4)),
    "ordinal"
  )
})

test_that("strictly decreasing scores do not produce a warning", {
  # regression test: decreasing scores are monotone and must not warn
  expect_no_warning(
    mantelTrendTest(Job, srow = c(4,3,2,1))
  )
})

test_that("decreasing scores give the same statistic as increasing (sign-invariant)", {
  res_inc <- suppressWarnings(mantelTrendTest(Job, srow = 1:4))
  res_dec <- suppressWarnings(mantelTrendTest(Job, srow = 4:1))

  expect_equal(unname(res_inc$statistic), unname(res_dec$statistic))
  expect_equal(res_inc$p.value, res_dec$p.value)
})
# -------------------------------------------------------------------------
# Automatic scores from numeric dimnames
# -------------------------------------------------------------------------
test_that("character dimnames (Job) still default to 1:nrow/1:ncol", {
  res_auto     <- mantelTrendTest(Job)
  res_explicit <- mantelTrendTest(Job, srow = 1:4, scol = 1:4)

  expect_equal(res_auto$statistic, res_explicit$statistic)
})

test_that("numeric row dimnames are used as default srow", {
  dose <- matrix(c(10,9,10,7, 0,1,0,3), nrow = 4,
                 dimnames = list(dose = c("0","1","5","20"),
                                 resp = c("no","yes")))

  res_auto   <- mantelTrendTest(dose)
  res_manual <- mantelTrendTest(dose, srow = c(0,1,5,20), scol = c(1,2))

  expect_equal(unname(res_auto$statistic), unname(res_manual$statistic))

  # and it must differ from the old 1:nrow(x) default, since these
  # dimnames are not evenly spaced
  res_old_default <- mantelTrendTest(dose, srow = 1:4, scol = 1:2)
  expect_false(isTRUE(all.equal(unname(res_auto$statistic),
                                unname(res_old_default$statistic))))
})

test_that("explicit srow/scol still override the automatic default", {
  dose <- matrix(c(10,9,10,7, 0,1,0,3), nrow = 4,
                 dimnames = list(dose = c("0","1","5","20"),
                                 resp = c("no","yes")))

  res <- mantelTrendTest(dose, srow = c(1, 2, 3, 4))

  expect_equal(unname(res$statistic),
               unname(mantelTrendTest(dose, srow = 1:4, scol = c(1,2))$statistic))
})

