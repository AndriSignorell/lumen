
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


test_that("exact inference now covers tied data", {
  
  set.seed(4)
  
  x <- sample(1:5, 15, TRUE)
  g <- rep(1:3, each = 5)
  
  res <- jonckheereTerpstraTest(x, g, method = "exact")
  
  expect_match(res$method, "exact, ties", fixed = TRUE)
  expect_silent(jonckheereTerpstraTest(x, g, method = "exact"))
  
  # the exact statistic is half-integral where tied pairs cross groups
  expect_equal(unname(res$statistic) %% 0.5, 0)
})


test_that("exact method falls back to asymptotic when ties get too costly", {
  
  # few ties in a larger sample: the table of group by distinct value has
  # more cells than the recursion is allowed to walk
  set.seed(3)
  
  x <- round(rnorm(60), 3)
  x[1:4] <- 0
  g <- rep(1:3, each = 20)
  
  expect_warning(
    res <- jonckheereTerpstraTest(x, g, method = "exact"),
    "falling back to the asymptotic approximation")
  
  expect_match(res$method, "asymptotic", fixed = TRUE)
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
  
  expect_equal(r_exact$p.value, r_perm$p.value, tolerance = 0.02)
})



# -------------------------------------------------------------------------
# Hollander & Wolfe reference example
# -------------------------------------------------------------------------
test_that("H&W example gives correct JT statistic and p-value", {
  
  x <- c(24, 61, 59, 98, 73, 68, 91, 94, 79, 63, 82, 89)
  g <- ordered(rep(1:3, c(3, 4, 5)))
  
  res <- jonckheereTerpstraTest(x, g, method = "exact", alternative="increasing")
  
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


# exact inference with ties; merge into test-jonckheereTerpstraTest.R

jtByHand <- function(x, g) {
  
  lev <- levels(g)
  J <- 0
  
  for (u in seq_len(nlevels(g) - 1L)) {
    for (v in (u + 1L):nlevels(g)) {
      a <- x[g == lev[u]]
      b <- x[g == lev[v]]
      J <- J + sum(outer(a, b, "<")) + 0.5 * sum(outer(a, b, "=="))
    }
  }
  
  J
}


# Hollander and Wolfe, Example 6.2
motiv <- list(
  no    = c(40, 35, 38, 43, 44, 41),
  rough = c(38, 40, 47, 44, 40, 42),
  acc   = c(48, 40, 45, 43, 46, 44))


test_that("the tied recursion reproduces the tie-free one", {
  
  gsize <- c(3L, 2L, 4L)
  
  tied <- .jtpdfTies(gsize, rep(1L, sum(gsize)))
  free <- .jtpdf(gsize)
  
  expect_equal(sum(tied), 1)
  
  # without ties the statistic is integral, the odd cells carry no mass
  expect_equal(tied[seq(1L, length(tied), by = 2L)], free)
  expect_true(all(tied[seq(2L, length(tied), by = 2L)] == 0))
})


test_that("the tied distribution is a distribution", {
  
  pdf <- .jtpdfTies(c(4L, 3L, 5L), c(2L, 3L, 2L, 1L, 4L))
  
  expect_equal(sum(pdf), 1)
  expect_true(all(pdf >= 0))
  
  # with equal group sizes reversing their order is a relabelling, so the
  # distribution is symmetric about its mean; unequal sizes are not
  equal <- .jtpdfTies(c(4L, 4L, 4L), c(2L, 3L, 2L, 1L, 4L))
  
  expect_equal(equal, rev(equal))
})


test_that("exact inference is used and reported for tied data", {
  
  res <- jonckheereTerpstraTest(motiv, alternative = "increasing")
  
  expect_equal(unname(res$statistic), 79)
  expect_equal(res$p.value, 0.02096649, tolerance = 1e-6)
  expect_match(res$method, "exact, ties", fixed = TRUE)
})


test_that("the reference values of Hollander and Wolfe are reproduced", {
  
  # Example 6.2, where the tie-free distribution and the tie-corrected
  # normal approximation are reported side by side, 0.0231 and 0.0207
  # rounded to the precision the book prints, the values being 0.02306371,
  # 0.02071039 and 150.2868
  pdf <- .jtpdf(c(6L, 6L, 6L))
  
  expect_equal(round(sum(pdf[80:length(pdf)]), 4), 0.0231)
  
  res <- jonckheereTerpstraTest(motiv, alternative = "increasing",
                                method = "asymptotic")
  
  expect_equal(round(res$p.value, 4), 0.0207)
  
  # the null variance of Equation 6.19, var_0(J) = 150.29
  z <- qnorm(res$p.value, lower.tail = FALSE)
  
  expect_equal(round((79 - 54)^2 / z^2, 2), 150.29)
})


test_that("the tie correction of the asymptotic variance is not enough", {
  
  # 15 of 16 observations tied, the single larger value in the last group
  x <- c(rep(57, 15), 59)
  g <- rep(1:4, 4)
  
  # the value can only land in one of the four groups, hence a quarter
  expect_equal(
    jonckheereTerpstraTest(x, g, alternative = "increasing")$p.value, 0.25)
  
  expect_equal(jonckheereTerpstraTest(x, g)$p.value, 0.5)
  
  # the normal approximation is far off in both directions here
  expect_equal(
    jonckheereTerpstraTest(x, g, method = "asymptotic")$p.value,
    0.1797, tolerance = 1e-4)
})


test_that("two groups reproduce the exact Wilcoxon rank sum test", {
  
  a <- c(1, 2, 2, 3, 5)
  b <- c(2, 4, 4, 6, 7)
  
  res <- jonckheereTerpstraTest(list(a = a, b = b), alternative = "increasing")
  
  # the exact conditional test of the tied data, which wilcox.test declines
  expect_match(res$method, "exact, ties", fixed = TRUE)
  
  untied <- jonckheereTerpstraTest(list(a = c(1, 2, 3), b = c(4, 5, 6)),
                                   alternative = "increasing")
  
  expect_equal(
    untied$p.value,
    wilcox.test(c(4, 5, 6), c(1, 2, 3), alternative = "greater")$p.value)
})


test_that("the exact p-value agrees with the permutation p-value", {
  
  set.seed(11)
  
  x <- sample(1:4, 13, TRUE)
  g <- rep(1:3, c(4, 4, 5))
  
  exact <- jonckheereTerpstraTest(x, g, alternative = "increasing")
  perm <- jonckheereTerpstraTest(x, g, alternative = "increasing",
                                 method = "permutation", R = 20000)
  
  expect_equal(exact$p.value, perm$p.value, tolerance = 0.01)
})


test_that("reversing the group order swaps the one-sided p-values", {
  
  x <- unlist(motiv, use.names = FALSE)
  g <- factor(rep(names(motiv), lengths(motiv)), levels = names(motiv))
  
  up <- jonckheereTerpstraTest(x, g, alternative = "increasing")
  down <- jonckheereTerpstraTest(x, factor(g, levels = rev(levels(g))),
                                 alternative = "decreasing")
  
  expect_equal(up$p.value, down$p.value)
})


test_that("constant data carry no evidence, whichever method is asked for", {
  
  x <- rep(5, 9)
  g <- rep(1:3, each = 3)
  
  res <- jonckheereTerpstraTest(x, g)
  
  expect_equal(unname(res$statistic), jtByHand(x, factor(g)))
  expect_equal(res$p.value, 1)
  
  # the asymptotic variance vanishes here, z would be 0/0
  expect_equal(jonckheereTerpstraTest(x, g, method = "asymptotic")$p.value, 1)
  expect_equal(
    jonckheereTerpstraTest(x, g, method = "permutation", R = 99)$p.value, 1)
  
  # and beyond the sample size at which the tie-free recursion stops
  big <- jonckheereTerpstraTest(rep(5, 150), rep(1:3, each = 50))
  
  expect_equal(big$p.value, 1)
  
  expect_equal(
    jonckheereTerpstraTest(rep(5, 150), rep(1:3, each = 50),
                           method = "asymptotic")$p.value, 1)
})


test_that("the three methods share one two-sided rule", {
  
  # smallest asymmetric case: JT is 0 with probability 1/3 and 1.5 with
  # probability 2/3, so the doubled lower tail is 2/3, while the values
  # farther from the null mean of 1 than the observed one are 1/3
  x <- c(1, 0, 0)
  g <- c(1, 2, 2)
  
  expect_equal(.jtpdfTies(c(1L, 2L), c(2L, 1L)), c(1/3, 0, 0, 2/3, 0))
  
  exact <- jonckheereTerpstraTest(x, g)
  
  expect_equal(exact$p.value, 2/3)
  
  set.seed(5)
  perm <- jonckheereTerpstraTest(x, g, method = "permutation", R = 20000)
  
  expect_equal(perm$p.value, 2/3, tolerance = 0.02)
  
  # one-sided the two agree anyway
  for (alt in c("increasing", "decreasing")) {
    set.seed(5)
    expect_equal(
      jonckheereTerpstraTest(x, g, alternative = alt)$p.value,
      jonckheereTerpstraTest(x, g, alternative = alt,
                             method = "permutation", R = 20000)$p.value,
      tolerance = 0.02)
  }
})


test_that("R is validated before any permutation is drawn", {
  
  x <- c(1, 0, 0)
  g <- c(1, 2, 2)
  
  for (bad in list(99.5, c(10, 20), 0, -5, NA_real_, Inf, "100")) {
    expect_error(
      jonckheereTerpstraTest(x, g, method = "permutation", R = bad),
      "single positive integer")
  }
  
  expect_error(jonckheereTerpstraTest(x, g, method = "permutation"),
               "must be specified")
})


test_that("the recursions reject impossible arguments", {
  
  expect_error(.jtpdfTies(c(2L, 2L), c(1L, 1L, 1L)), "sum to the number")
  expect_error(.jtpdfTies(c(3L, 0L), c(1L, 1L, 1L)), "positive")
  expect_error(.jtpdf(c(3L, 0L)), "positive")
  
  # a table that could not be held in memory is refused before the
  # recursion is entered at all
  expect_error(.jtpdfTies(c(20L, 20L), rep(1L, 40)), "too large")
  expect_gt(.jtTiesCells(c(20L, 20L), rep(1L, 40)), .jtTiesMaxCells)
})


test_that("the cost of the tied recursion decides which method is used", {
  
  # few ties in a large sample: the table has more cells than can be walked
  set.seed(3)
  
  x <- round(rnorm(60), 3)
  x[1:4] <- 0
  g <- rep(1:3, each = 20)
  
  expect_gt(.jtTiesCells(rep(20L, 3), as.integer(table(x))), .jtTiesMaxCells)
  
  expect_warning(
    res <- jonckheereTerpstraTest(x, g, method = "exact"),
    "falling back to the asymptotic approximation")
  
  expect_match(res$method, "asymptotic", fixed = TRUE)
  
  # auto turns to the approximation well before that limit
  expect_match(jonckheereTerpstraTest(x, g)$method, "asymptotic", fixed = TRUE)
  
  # the same sample size with three distinct values only stays within reach
  expect_lt(.jtTiesCells(rep(20L, 3), c(20L, 20L, 20L)), .jtTiesMaxCells)
})


test_that(".jtTiesCells counts states times support", {
  
  gsize <- c(6L, 6L, 6L)
  cnt <- as.integer(table(unlist(motiv)))
  
  maxJ <- 6 * 12 + 6 * 6
  
  expect_equal(.jtTiesCells(gsize, cnt), prod(cnt + 1) * (2 * maxJ + 1))
  expect_length(.jtpdfTies(gsize, cnt), 2 * maxJ + 1)
})
