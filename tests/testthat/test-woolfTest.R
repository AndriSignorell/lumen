

migraine <- xtabs(freq ~ .,
                  cbind(expand.grid(treatment = c("active","placebo"),
                                    response  = c("better","same"),
                                    gender    = c("female","male")),
                        freq = c(16,5,11,20,12,7,16,19)))
# -------------------------------------------------------------------------
# Basic functionality
# -------------------------------------------------------------------------
test_that("woolfTest returns an htest object", {
  
  res <- woolfTest(migraine)
  
  expect_s3_class(res, "htest")
  expect_named(res$statistic, "X-squared")
  expect_named(res$parameter, "df")
  expect_true(is.numeric(res$p.value))
  expect_gte(res$p.value, 0)
  expect_lte(res$p.value, 1)
})
# -------------------------------------------------------------------------
# Reference value (vcd)
# -------------------------------------------------------------------------
test_that("migraine example matches vcd::woolf_test", {
  
  res <- woolfTest(migraine)
  
  expect_equal(unname(res$statistic["X-squared"]), 1.4808, tolerance = 1e-3)
  expect_equal(res$p.value, 0.2236, tolerance = 1e-3)
  expect_equal(unname(res$parameter["df"]), 1L)
})
# -------------------------------------------------------------------------
# Zero cell correction
# -------------------------------------------------------------------------
test_that("zero cell triggers message and adds 0.5", {
  
  x <- migraine
  x[1,1,1] <- 0
  
  expect_message(woolfTest(x), "0.5")
})

test_that("zero cell result differs from no-zero result", {
  
  res_clean <- woolfTest(migraine)
  
  x_zero        <- migraine
  x_zero[1,1,1] <- 0
  
  res_zero <- suppressMessages(woolfTest(x_zero))
  
  expect_false(
    isTRUE(all.equal(
      unname(res_clean$statistic),
      unname(res_zero$statistic)
    ))
  )
})
# -------------------------------------------------------------------------
# observed and expected in output
# -------------------------------------------------------------------------
test_that("observed log-ORs have length K", {
  
  res <- woolfTest(migraine)
  
  expect_length(res$observed, dim(migraine)[3L])
})

test_that("expected is a scalar weighted mean", {
  
  res <- woolfTest(migraine)
  
  expect_length(res$expected, 1L)
  expect_true(is.numeric(res$expected))
})
# -------------------------------------------------------------------------
# df = K - 1
# -------------------------------------------------------------------------
test_that("df equals K - 1 for 3 strata", {
  
  x3 <- array(c(10,5,8,12, 8,3,6,9, 12,4,10,15),
              dim = c(2,2,3))
  
  res <- woolfTest(x3)
  
  expect_equal(unname(res$parameter["df"]), 2L)
})
# -------------------------------------------------------------------------
# Homogeneous ORs give non-significant result
# -------------------------------------------------------------------------
test_that("identical strata give X-squared near 0", {
  
  stratum <- matrix(c(20,10,5,15), nrow = 2)
  x_hom   <- array(c(stratum, stratum), dim = c(2,2,2))
  
  res <- woolfTest(x_hom)
  
  expect_lt(unname(res$statistic), 1e-10)
  expect_gt(res$p.value, 0.99)
})
# -------------------------------------------------------------------------
# Input validation
# -------------------------------------------------------------------------
test_that("non-array input throws error", {
  expect_error(woolfTest(matrix(1:4, 2)), "2x2xK")
})

test_that("non-2x2 strata throw error", {
  expect_error(woolfTest(array(1:27, dim = c(3,3,3))), "2x2xK")
})

test_that("negative counts throw error", {
  x_neg        <- migraine
  x_neg[1,1,1] <- -1
  expect_error(woolfTest(x_neg), "nonnegative")
})

test_that("non-integer counts produce warning", {
  x_frac        <- migraine
  x_frac[1,1,1] <- 16.5
  expect_warning(woolfTest(x_frac), "non-integer")
})
# -------------------------------------------------------------------------
# Consistent with breslowDayTest direction
# -------------------------------------------------------------------------
test_that("woolfTest and breslowDayTest agree on significance direction", {
  
  res_w <- woolfTest(migraine)
  res_b <- breslowDayTest(migraine)
  
  # Both should give non-significant result for migraine
  expect_gt(res_w$p.value, 0.05)
  expect_gt(res_b$p.value, 0.05)
})
# -------------------------------------------------------------------------
# Print compatibility
# -------------------------------------------------------------------------
test_that("print.htest works", {
  
  res <- woolfTest(migraine)
  
  expect_output(print(res), "Woolf")
})

