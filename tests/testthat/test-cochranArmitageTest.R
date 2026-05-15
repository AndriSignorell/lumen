
# -------------------------------------------------------------------------
# Basic functionality
# -------------------------------------------------------------------------
test_that("cochranArmitageTest returns an htest object", {
  
  dose <- matrix(c(10,9,10,7, 0,1,0,3), byrow = TRUE, nrow = 2,
                 dimnames = list(resp = 0:1, dose = 0:3))
  
  res <- cochranArmitageTest(dose)
  
  expect_s3_class(res, "htest")
  expect_named(res$statistic, "Z")
  expect_named(res$parameter, "dim")
  expect_true(is.numeric(res$p.value))
  expect_gte(res$p.value, 0)
  expect_lte(res$p.value, 1)
})
# -------------------------------------------------------------------------
# Known reference values
# -------------------------------------------------------------------------
test_that("lungtumor example gives correct Z and p-value", {
  
  lungtumor <- data.frame(
    dose  = rep(c(0,1,2), c(40,50,48)),
    tumor = c(rep(c(0,1), c(38,2)),
              rep(c(0,1), c(43,7)),
              rep(c(0,1), c(33,15)))
  )
  
  tab <- table(lungtumor$dose, lungtumor$tumor)
  res <- cochranArmitageTest(tab)
  
  expect_equal(unname(res$statistic["Z"]), -3.2735, tolerance = 1e-5)
  expect_equal(res$p.value, 0.001062295, tolerance = 1e-6)
})

test_that("pain example gives correct Z and p-value (SAS reference)", {
  
  pain <- structure(c(26,6,26,7,23,9,18,14,9,23),
                    dim = c(2L,5L),
                    dimnames = list(adverse = c("No","Yes"),
                                    dose    = c("0","1","2","3","4")),
                    class = "table")
  
  res <- cochranArmitageTest(pain)
  
  # SAS PROC FREQ Output 41.8.4: Z = -4.7918, p < .0001
  expect_equal(unname(res$statistic["Z"]), -4.7918, tolerance = 1e-3)
  expect_lt(res$p.value, 0.0001)
})
# -------------------------------------------------------------------------
# Transposition invariance
# -------------------------------------------------------------------------
test_that("transposed input gives same result", {
  
  pain <- structure(c(26,6,26,7,23,9,18,14,9,23),
                    dim = c(2L,5L),
                    dimnames = list(adverse = c("No","Yes"),
                                    dose    = c("0","1","2","3","4")),
                    class = "table")
  
  res1 <- cochranArmitageTest(pain)
  res2 <- cochranArmitageTest(t(pain))
  
  expect_equal(unname(res1$statistic), unname(res2$statistic))
  expect_equal(res1$p.value, res2$p.value)
})
# -------------------------------------------------------------------------
# Alternatives
# -------------------------------------------------------------------------
test_that("decreasing alternative detects negative trend in pain example", {
  
  pain <- structure(c(26,6,26,7,23,9,18,14,9,23),
                    dim = c(2L,5L),
                    dimnames = list(adverse = c("No","Yes"),
                                    dose    = c("0","1","2","3","4")),
                    class = "table")
  
  res <- cochranArmitageTest(pain, alternative = "decreasing")
  
  expect_lt(res$p.value, 0.05)
})

test_that("increasing alternative detects positive trend in lungtumor example", {
  
  lungtumor <- data.frame(
    dose  = rep(c(0,1,2), c(40,50,48)),
    tumor = c(rep(c(0,1), c(38,2)),
              rep(c(0,1), c(43,7)),
              rep(c(0,1), c(33,15)))
  )
  
  # tumor = 0 (no tumor) decreases with dose → "decreasing"
  tab <- table(lungtumor$dose, lungtumor$tumor)
  res <- cochranArmitageTest(tab, alternative = "decreasing")
  
  expect_lt(res$p.value, 0.05)
})

test_that("two-sided p-value >= one-sided p-value", {
  
  pain <- structure(c(26,6,26,7,23,9,18,14,9,23),
                    dim = c(2L,5L),
                    dimnames = list(adverse = c("No","Yes"),
                                    dose    = c("0","1","2","3","4")),
                    class = "table")
  
  r_dec <- cochranArmitageTest(pain, alternative = "decreasing")
  r_two <- cochranArmitageTest(pain, alternative = "two.sided")
  
  expect_gte(r_two$p.value, r_dec$p.value)
})

test_that("increasing and decreasing p-values sum to 1", {
  
  pain <- structure(c(26,6,26,7,23,9,18,14,9,23),
                    dim = c(2L,5L),
                    dimnames = list(adverse = c("No","Yes"),
                                    dose    = c("0","1","2","3","4")),
                    class = "table")
  
  r_inc <- cochranArmitageTest(pain, alternative = "increasing")
  r_dec <- cochranArmitageTest(pain, alternative = "decreasing")
  
  expect_equal(r_inc$p.value + r_dec$p.value, 1, tolerance = 1e-10)
})

test_that("two-sided p = 2 * min(increasing, decreasing)", {
  
  pain <- structure(c(26,6,26,7,23,9,18,14,9,23),
                    dim = c(2L,5L),
                    dimnames = list(adverse = c("No","Yes"),
                                    dose    = c("0","1","2","3","4")),
                    class = "table")
  
  r_inc <- cochranArmitageTest(pain, alternative = "increasing")
  r_dec <- cochranArmitageTest(pain, alternative = "decreasing")
  r_two <- cochranArmitageTest(pain, alternative = "two.sided")
  
  expect_equal(
    r_two$p.value,
    2 * min(r_inc$p.value, r_dec$p.value),
    tolerance = 1e-10
  )
})

test_that("invalid alternative throws error", {
  
  tab <- matrix(c(10,9,10,7, 0,1,0,3), byrow = TRUE, nrow = 2)
  
  expect_error(
    cochranArmitageTest(tab, alternative = "invalid")
  )
})
# -------------------------------------------------------------------------
# Numeric dimnames used as scores
# -------------------------------------------------------------------------
test_that("numeric dimnames produce same Z as sequential scores when equidistant", {
  
  # Scores 0,1,2,3 vs 1,2,3,4: nur Verschiebung, Rbar passt sich an → Z identisch
  tab_named   <- matrix(c(10,9,10,7, 0,1,0,3), byrow = TRUE, nrow = 2,
                        dimnames = list(resp = 0:1, dose = 0:3))
  tab_unnamed <- matrix(c(10,9,10,7, 0,1,0,3), byrow = TRUE, nrow = 2)
  
  res_named   <- cochranArmitageTest(tab_named)
  res_unnamed <- cochranArmitageTest(tab_unnamed)
  
  expect_equal(
    unname(res_named$statistic),
    unname(res_unnamed$statistic),
    tolerance = 1e-10
  )
})
# -------------------------------------------------------------------------
# Input validation
# -------------------------------------------------------------------------
test_that("non-2-level table throws error", {
  
  tab <- matrix(c(10,9,10, 5,6,7, 1,2,3), nrow = 3)
  
  expect_error(
    cochranArmitageTest(tab),
    "rx2"
  )
})

test_that("parameter dim reflects number of ordered levels", {
  
  tab <- matrix(c(10,9,10,7, 0,1,0,3), byrow = TRUE, nrow = 2)
  
  res <- cochranArmitageTest(tab)
  
  expect_equal(unname(res$parameter["dim"]), 4L)
})
# -------------------------------------------------------------------------
# Print compatibility
# -------------------------------------------------------------------------
test_that("print.htest works", {
  
  tab <- matrix(c(10,9,10,7, 0,1,0,3), byrow = TRUE, nrow = 2)
  
  res <- cochranArmitageTest(tab)
  
  expect_output(print(res), "Cochran-Armitage")
})

