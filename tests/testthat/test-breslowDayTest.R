

# -------------------------------------------------------------------------
# Setup
# -------------------------------------------------------------------------
migraine <- xtabs(freq ~ .,
                  cbind(expand.grid(treatment = c("active","placebo"),
                                    response  = c("better","same"),
                                    gender    = c("female","male")),
                        freq = c(16,5,11,20,12,7,16,19)))

salary <- array(
  c(38,12,102,141,12,9,136,383),
  dim = c(2,2,2),
  dimnames = list(
    exposure = c("exposed","not"),
    disease  = c("case","control"),
    salary   = c("<1000",">=1000")
  )
)
# -------------------------------------------------------------------------
# Basic functionality
# -------------------------------------------------------------------------
test_that("breslowDayTest returns an htest object", {
  
  res <- breslowDayTest(migraine)
  
  expect_s3_class(res, "htest")
  expect_named(res$statistic, "X-squared")
  expect_named(res$parameter, "df")
  expect_true(is.numeric(res$p.value))
  expect_gte(res$p.value, 0)
  expect_lte(res$p.value, 1)
})
# -------------------------------------------------------------------------
# Reference values (migraine)
# -------------------------------------------------------------------------
test_that("migraine example gives correct statistic and p-value", {
  
  res <- breslowDayTest(migraine)
  
  expect_equal(unname(res$statistic["X-squared"]), 1.4929, tolerance = 1e-3)
  expect_equal(res$p.value, 0.2218, tolerance = 1e-3)
  expect_equal(unname(res$parameter["df"]), 1L)
})

test_that("migraine Tarone correction gives correct statistic", {
  
  res <- breslowDayTest(migraine, correct = TRUE)
  
  expect_equal(unname(res$statistic["X-squared"]), 1.4905, tolerance = 1e-3)
  expect_equal(res$p.value, 0.2221, tolerance = 1e-3)
})
# -------------------------------------------------------------------------
# Tarone correction reduces statistic
# -------------------------------------------------------------------------
test_that("Tarone correction gives smaller or equal statistic", {
  
  res_plain   <- breslowDayTest(migraine)
  res_tarone  <- breslowDayTest(migraine, correct = TRUE)
  
  expect_lte(
    unname(res_tarone$statistic),
    unname(res_plain$statistic)
  )
})
# -------------------------------------------------------------------------
# Custom OR
# -------------------------------------------------------------------------
test_that("custom OR gives different result than MH estimate", {
  
  res_mh  <- breslowDayTest(salary)
  res_or  <- breslowDayTest(salary, OR = 4.02)
  
  expect_false(
    isTRUE(all.equal(
      unname(res_mh$statistic),
      unname(res_or$statistic)
    ))
  )
})

test_that("OR = MH estimate gives same result as default", {
  
  # Compute MH estimate manually
  a <- salary[1,1,]; b <- salary[1,2,]
  c <- salary[2,1,]; d <- salary[2,2,]
  n <- a + b + c + d
  or_mh <- sum(a*d/n) / sum(b*c/n)
  
  res_default <- breslowDayTest(salary)
  res_manual  <- breslowDayTest(salary, OR = or_mh)
  
  expect_equal(
    unname(res_default$statistic),
    unname(res_manual$statistic),
    tolerance = 1e-10
  )
})
# -------------------------------------------------------------------------
# df = K - 1
# -------------------------------------------------------------------------
test_that("df equals K - 1", {
  
  # 3 strata
  x3 <- array(c(10,5,8,12, 8,3,6,9, 12,4,10,15),
              dim = c(2,2,3))
  
  res <- breslowDayTest(x3)
  
  expect_equal(unname(res$parameter["df"]), 2L)
})
# -------------------------------------------------------------------------
# Homogeneous ORs give non-significant result
# -------------------------------------------------------------------------
test_that("identical strata give X-squared near 0", {
  
  # Both strata identical -> ORs identical -> no heterogeneity
  stratum <- matrix(c(20,10,5,15), nrow = 2)
  x_hom   <- array(c(stratum, stratum), dim = c(2,2,2))
  
  res <- breslowDayTest(x_hom)
  
  expect_lt(unname(res$statistic), 1e-6)
  expect_gt(res$p.value, 0.99)
})
# -------------------------------------------------------------------------
# Input validation
# -------------------------------------------------------------------------
test_that("non-array input throws error", {
  
  expect_error(breslowDayTest(matrix(1:4, 2)), "2x2xK")
})

test_that("non-2x2 strata throw error", {
  
  x <- array(1:27, dim = c(3,3,3))
  
  expect_error(breslowDayTest(x), "2x2xK")
})

test_that("method string reflects correct argument", {
  
  res_plain  <- breslowDayTest(migraine)
  res_tarone <- breslowDayTest(migraine, correct = TRUE)
  
  expect_false(grepl("Tarone", res_plain$method))
  expect_match(res_tarone$method, "Tarone")
})
# -------------------------------------------------------------------------
# statistic name is clean
# -------------------------------------------------------------------------
test_that("statistic name is X-squared without suffix", {
  
  res <- breslowDayTest(migraine)
  
  expect_equal(names(res$statistic), "X-squared")
})
# -------------------------------------------------------------------------
# Print compatibility
# -------------------------------------------------------------------------
test_that("print.htest works", {
  
  res <- breslowDayTest(migraine)
  
  expect_output(print(res), "Breslow")
})


# -------------------------------------------------------------------------
# Input validation – neue Tests
# -------------------------------------------------------------------------
test_that("invalid OR throws error", {
  expect_error(breslowDayTest(migraine, OR = -1),   "positive")
  expect_error(breslowDayTest(migraine, OR = 0),    "positive")
  expect_error(breslowDayTest(migraine, OR = Inf),  "positive")
  expect_error(breslowDayTest(migraine, OR = "x"),  "positive")
})

test_that("negative counts throw error", {
  x_neg <- migraine
  x_neg[1,1,1] <- -1
  expect_error(breslowDayTest(x_neg), "nonnegative")
})

test_that("non-integer counts produce warning", {
  x_frac <- migraine
  x_frac[1,1,1] <- 16.5
  expect_warning(breslowDayTest(x_frac), "non-integer")
})

test_that("invalid correct throws error", {
  expect_error(breslowDayTest(migraine, correct = NA), "TRUE or FALSE")
})

test_that("n equals total count", {
  res <- breslowDayTest(migraine)
  expect_equal(res$n, sum(migraine))
})

