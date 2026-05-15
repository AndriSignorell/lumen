
# -------------------------------------------------------------------------
# Basic functionality
# -------------------------------------------------------------------------

test_that("jonckheereTerpstraTest returns an htest object", {
  
  set.seed(1)
  
  g <- ordered(rep(1:4, each = 10))
  x <- rnorm(40)
  
  res <- jonckheereTerpstraTest(x, g)
  
  expect_s3_class(res, "htest")
  expect_named(res$parameter, c("k", "n"))
  expect_named(res$statistic, "JT")
  expect_true(is.numeric(res$p.value))
  expect_true(res$p.value >= 0)
  expect_true(res$p.value <= 1)
})
# -------------------------------------------------------------------------
# Increasing trend detection
# -------------------------------------------------------------------------
test_that("detects increasing trend", {
  
  set.seed(1)
  
  g <- ordered(rep(1:5, each = 20))
  x <- rnorm(100) + as.numeric(g)
  
  res <- jonckheereTerpstraTest(x, g, alternative = "increasing")
  
  expect_lt(res$p.value, 0.05)
})
# -------------------------------------------------------------------------
# Decreasing trend detection
# -------------------------------------------------------------------------
test_that("detects decreasing trend", {
  
  set.seed(1)
  
  g <- ordered(rep(1:5, each = 20))
  x <- rnorm(100) - as.numeric(g)
  
  res <- jonckheereTerpstraTest(x, g, alternative = "decreasing")
  
  expect_lt(res$p.value, 0.05)
})
# -------------------------------------------------------------------------
# No trend: p-value not systematically small
# -------------------------------------------------------------------------
test_that("no trend gives non-significant p-value", {
  
  set.seed(42)
  
  g <- ordered(rep(1:4, each = 25))
  x <- rnorm(100)
  
  res <- jonckheereTerpstraTest(x, g, alternative = "increasing")
  
  expect_gt(res$p.value, 0.05)
})
# -------------------------------------------------------------------------
# Two-sided alternative
# -------------------------------------------------------------------------
test_that("two-sided alternative works", {
  
  set.seed(1)
  
  g <- ordered(rep(1:4, each = 15))
  x <- rnorm(60)
  
  res <- jonckheereTerpstraTest(x, g, alternative = "two.sided")
  
  expect_true(is.finite(res$p.value))
  expect_gte(res$p.value, 0)
  expect_lte(res$p.value, 1)
})

test_that("two-sided p-value >= one-sided p-value", {
  
  set.seed(1)
  
  g <- ordered(rep(1:4, each = 10))
  x <- rnorm(40) + as.numeric(g)
  
  r_inc <- jonckheereTerpstraTest(x, g, alternative = "increasing")
  r_two <- jonckheereTerpstraTest(x, g, alternative = "two.sided")
  
  expect_gte(r_two$p.value, r_inc$p.value)
})
# -------------------------------------------------------------------------
# Formula interface
# -------------------------------------------------------------------------
test_that("formula interface gives same result as default", {
  
  set.seed(1)
  
  g   <- ordered(rep(1:4, each = 10))
  x   <- rnorm(40)
  dat <- data.frame(x = x, g = g)
  
  res1 <- jonckheereTerpstraTest(x, g)
  res2 <- jonckheereTerpstraTest(x ~ g, data = dat)
  
  expect_equal(unname(res1$statistic), unname(res2$statistic))
  expect_equal(unname(res1$p.value),   unname(res2$p.value))
})
# -------------------------------------------------------------------------
# List interface
# -------------------------------------------------------------------------
test_that("list interface works", {
  
  set.seed(1)
  
  x <- list(
    a = rnorm(10),
    b = rnorm(10, 1),
    c = rnorm(10, 2)
  )
  
  res <- jonckheereTerpstraTest(x)
  
  expect_s3_class(res, "htest")
  expect_equal(unname(res$parameter["k"]), 3)
  expect_equal(unname(res$parameter["n"]), 30)
})

test_that("list interface gives same result as default interface", {
  
  set.seed(1)
  
  xa <- rnorm(8)
  xb <- rnorm(8, 1)
  xc <- rnorm(8, 2)
  
  x <- c(xa, xb, xc)
  g <- ordered(rep(1:3, each = 8))
  
  res_default <- jonckheereTerpstraTest(x, g)
  res_list    <- jonckheereTerpstraTest(list(xa, xb, xc))
  
  expect_equal(
    unname(res_default$statistic),
    unname(res_list$statistic)
  )
})
# -------------------------------------------------------------------------
# method = "exact"
# -------------------------------------------------------------------------
test_that("exact method works without ties", {
  
  set.seed(1)
  
  g <- ordered(rep(1:3, each = 5))
  x <- rnorm(15)
  
  res <- jonckheereTerpstraTest(x, g, method = "exact")
  
  expect_match(res$method, "exact")
  expect_true(is.finite(res$p.value))
  expect_gte(res$p.value, 0)
  expect_lte(res$p.value, 1)
})

test_that("exact method falls back to asymptotic with ties", {
  
  g <- ordered(rep(1:3, each = 5))
  x <- c(1, 1, 2, 3, 4,
         2, 2, 3, 4, 5,
         3, 3, 4, 5, 6)
  
  expect_warning(
    res <- jonckheereTerpstraTest(x, g, method = "exact"),
    "ties"
  )
  
  expect_match(res$method, "asymptotic")
})

test_that("exact method warns for n > 100", {
  
  g <- ordered(rep(1:3, each = 40))
  x <- seq_len(120) + rnorm(120)
  
  expect_warning(
    jonckheereTerpstraTest(x, g, method = "exact"),
    "100"
  )
})

test_that("auto method selects exact for small n without ties", {
  
  set.seed(1)
  
  g <- ordered(rep(1:3, each = 5))
  x <- rnorm(15)
  
  res <- jonckheereTerpstraTest(x, g, method = "auto")
  
  expect_match(res$method, "exact")
})

test_that("auto method selects asymptotic for large n", {
  
  set.seed(1)
  
  g <- ordered(rep(1:3, each = 40))
  x <- rnorm(120)
  
  res <- jonckheereTerpstraTest(x, g, method = "auto")
  
  expect_match(res$method, "asymptotic")
})
# -------------------------------------------------------------------------
# method = "asymptotic"
# -------------------------------------------------------------------------
test_that("asymptotic method works", {
  
  set.seed(1)
  
  g <- ordered(rep(1:4, each = 10))
  x <- rnorm(40)
  
  res <- jonckheereTerpstraTest(x, g, method = "asymptotic")
  
  expect_match(res$method, "asymptotic")
  expect_true(is.finite(res$p.value))
})
# -------------------------------------------------------------------------
# method = "permutation"
# -------------------------------------------------------------------------
test_that("permutation method works", {
  
  set.seed(1)
  
  g <- ordered(rep(1:4, each = 10))
  x <- rnorm(40)
  
  res <- jonckheereTerpstraTest(x, g, method = "permutation", R = 500)
  
  expect_match(res$method, "permutation")
  expect_match(res$method, "500")
  expect_gte(res$p.value, 0)
  expect_lte(res$p.value, 1)
})

test_that("permutation method requires R", {
  
  g <- ordered(rep(1:3, each = 10))
  x <- rnorm(30)
  
  expect_error(
    jonckheereTerpstraTest(x, g, method = "permutation"),
    "R"
  )
})

test_that("R is ignored with warning when method != permutation", {
  
  set.seed(1)
  
  g <- ordered(rep(1:3, each = 10))
  x <- rnorm(30)
  
  expect_warning(
    jonckheereTerpstraTest(x, g, method = "asymptotic", R = 500),
    "ignored"
  )
})

test_that("permutation p-values stay within [0,1] across replications", {
  
  set.seed(1)
  
  for (i in seq_len(20)) {
    
    g <- ordered(rep(1:4, each = 10))
    x <- rnorm(40)
    
    res <- jonckheereTerpstraTest(x, g, method = "permutation", R = 200)
    
    expect_true(is.finite(res$p.value))
    expect_gte(res$p.value, 0)
    expect_lte(res$p.value, 1)
  }
})
# -------------------------------------------------------------------------
# Exact and asymptotic agree directionally
# -------------------------------------------------------------------------
test_that("exact and asymptotic agree directionally on strong trend", {
  
  set.seed(1)
  
  g <- ordered(rep(1:4, each = 8))
  x <- rnorm(32) + as.numeric(g)
  
  r_exact <- jonckheereTerpstraTest(x, g, method = "exact")
  r_asymp <- jonckheereTerpstraTest(x, g, method = "asymptotic")
  
  expect_lt(r_exact$p.value, 0.05)
  expect_lt(r_asymp$p.value, 0.05)
})


test_that("exact and permutation p-values are close", {

  set.seed(1)
    
  g <- ordered(rep(1:3, c(3, 4, 5)))
  x <- c(24, 61, 59, 98, 73, 68, 91, 94, 79, 63, 82, 89)
  
  r_exact <- jonckheereTerpstraTest(x, g, method = "exact")
  
  set.seed(1)
  r_perm <- jonckheereTerpstraTest(x, g, method = "permutation", R = 100000)
  
  expect_equal(r_exact$p.value, r_perm$p.value, tolerance = 0.01)
})



# -------------------------------------------------------------------------
# Hollander & Wolfe reference example
# -------------------------------------------------------------------------
test_that("H&W example gives correct JT statistic and p-value", {
  
  x <- c(24, 61, 59, 98, 73, 68, 91, 94, 79, 63, 82, 89)
  g <- ordered(rep(1:3, c(3, 4, 5)))
  
  res <- jonckheereTerpstraTest(x, g, method = "exact")
  
  expect_equal(unname(res$statistic["JT"]), 36)
  expect_equal(res$p.value, 0.03791, tolerance = 1e-3)
})


# -------------------------------------------------------------------------
# Missing values and non-finite values
# -------------------------------------------------------------------------
test_that("missing values are removed", {
  
  g <- ordered(rep(1:3, each = 5))
  x <- c(rnorm(13), NA, Inf)
  
  res <- jonckheereTerpstraTest(x, g)
  
  expect_true(is.finite(res$p.value))
})

test_that("result is unchanged after manually removing NAs", {
  
  set.seed(1)
  
  g_full <- ordered(rep(1:3, each = 6))
  x_full <- rnorm(18)
  x_full[c(2, 11)] <- NA
  
  ok      <- !is.na(x_full)
  res_na  <- jonckheereTerpstraTest(x_full, g_full)
  res_clean <- jonckheereTerpstraTest(x_full[ok], g_full[ok])
  
  expect_equal(
    unname(res_na$statistic),
    unname(res_clean$statistic)
  )
})
# -------------------------------------------------------------------------
# Input validation
# -------------------------------------------------------------------------
test_that("group length mismatch throws error", {
  expect_error(
    jonckheereTerpstraTest(1:10, factor(1:5)),
    "same length"
  )
})

test_that("single-group data throws error", {
  expect_error(
    jonckheereTerpstraTest(1:10, factor(rep(1, 10))),
    "same group"
  )
})

test_that("too few observations throws error", {
  expect_error(
    jonckheereTerpstraTest(numeric(1), factor(1)),
    "not enough"
  )
})

test_that("invalid alternative throws error", {
  expect_error(
    jonckheereTerpstraTest(rnorm(20), ordered(rep(1:2, each = 10)),
                           alternative = "invalid")
  )
})

test_that("invalid method throws error", {
  expect_error(
    jonckheereTerpstraTest(rnorm(20), ordered(rep(1:2, each = 10)),
                           method = "invalid")
  )
})
# -------------------------------------------------------------------------
# Invalid list input
# -------------------------------------------------------------------------
test_that("list with one group throws error", {
  expect_error(
    jonckheereTerpstraTest(list(1:5)),
    "at least two"
  )
})

test_that("list with empty group throws error", {
  expect_error(
    jonckheereTerpstraTest(list(a = numeric(0), b = 1:5)),
    "contain observations"
  )
})
# -------------------------------------------------------------------------
# Print compatibility
# -------------------------------------------------------------------------
test_that("print.htest works", {
  
  set.seed(1)
  
  g <- ordered(rep(1:3, each = 10))
  x <- rnorm(30)
  
  res <- jonckheereTerpstraTest(x, g)
  
  expect_output(print(res), "Jonckheere")
})

