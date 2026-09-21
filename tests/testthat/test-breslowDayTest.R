

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



test_that("OR = NA is backward compatible with OR = NULL (MH estimate)", {

  res_null <- breslowDayTest(salary)
  res_na   <- breslowDayTest(salary, OR = NA)

  expect_equal(res_null$statistic, res_na$statistic, tolerance = 1e-12)
  expect_equal(res_null$p.value, res_na$p.value, tolerance = 1e-12)
})


test_that("Tarone correction with user-supplied OR warns", {

  expect_warning(
    breslowDayTest(salary, OR = 4.02, correct = TRUE),
    "Tarone"
  )
})


# -- added --------------------------------------------------------------------

# expected count of cell (1,1) under a common OR, found by root search
tildeA <- function(tab, or) {
  m1 <- sum(tab[1, ]); n1 <- sum(tab[, 1]); N <- sum(tab)
  lo <- max(0, m1 + n1 - N); hi <- min(m1, n1)
  f <- function(a) a * (N - m1 - n1 + a) - or * (m1 - a) * (n1 - a)
  uniroot(f, c(lo, hi), tol = 1e-12)$root
}

bdRef <- function(x, or) {
  s <- 0; a <- ta <- va <- numeric(dim(x)[3])
  for (j in seq_len(dim(x)[3])) {
    t <- x[, , j]; e <- tildeA(t, or)
    m1 <- sum(t[1, ]); n1 <- sum(t[, 1]); N <- sum(t)
    v <- 1 / (1 / e + 1 / (m1 - e) + 1 / (n1 - e) + 1 / (N - m1 - n1 + e))
    s <- s + (t[1, 1] - e)^2 / v
    a[j] <- t[1, 1]; ta[j] <- e; va[j] <- v
  }
  list(stat = s, tarone = s - (sum(a) - sum(ta))^2 / sum(va))
}

orMH <- function(x) {
  n <- apply(x, 3, sum)
  sum(x[1, 1, ] * x[2, 2, ] / n) / sum(x[1, 2, ] * x[2, 1, ] / n)
}

test_that("statistic equals the definition with a root-searched expectation", {
  for (x in list(migraine, salary)) {
    ref <- bdRef(x, orMH(x))
    expect_equal(unname(breslowDayTest(x)$statistic), ref$stat, tolerance = 1e-8)
    expect_equal(unname(breslowDayTest(x, correct = TRUE)$statistic),
                 ref$tarone, tolerance = 1e-8)
  }
})

test_that("OR = 1 takes the linear branch: expectation m1 * n1 / N", {
  ref <- bdRef(salary, 1)
  expect_equal(unname(breslowDayTest(salary, OR = 1)$statistic), ref$stat,
               tolerance = 1e-8)
})

test_that("a hypothesised OR is not estimated: K df instead of K - 1", {
  r <- breslowDayTest(salary, OR = 4.02)
  expect_equal(unname(r$parameter), 2L)
  expect_equal(r$p.value, pchisq(unname(r$statistic), 2, lower.tail = FALSE))
  # OR = NA is the MH estimate, i.e. still K - 1
  expect_equal(unname(breslowDayTest(salary, OR = NA)$parameter), 1L)
})

test_that("size under H0 with a hypothesised OR", {
  set.seed(1)
  K <- 4
  p <- replicate(1500, {
    x <- array(0, c(2, 2, K))
    for (j in seq_len(K)) {
      p0 <- 0.3; p1 <- 2 * p0 / (1 - p0 + 2 * p0)
      x[, 1, j] <- c(rbinom(1, 60, p1), rbinom(1, 60, p0))
      x[, 2, j] <- 60 - x[, 1, j]
    }
    breslowDayTest(x, OR = 2)$p.value
  })
  expect_lt(abs(mean(p < 0.05) - 0.05), 0.025)
})

test_that("input checks", {
  x <- salary
  expect_error(breslowDayTest(array(1:8, c(2, 2, 2, 1))), "2x2xK")
  x[1, 1, 1] <- NA
  expect_error(breslowDayTest(x), "nonnegative and finite")
  expect_error(breslowDayTest(salary, OR = c(1, 2)), "positive finite")
  expect_error(breslowDayTest(salary, correct = c(TRUE, FALSE)), "TRUE or FALSE")
  expect_error(breslowDayTest(salary, correct = "a"), "TRUE or FALSE")

  z <- salary; z[1, , 1] <- 0
  expect_error(breslowDayTest(z), "zero marginal totals")

  # all b*c products zero: MH estimate undefined
  y <- array(c(5, 0, 0, 7, 4, 0, 0, 9), c(2, 2, 2))
  expect_error(breslowDayTest(y), "denominator is zero")
})
