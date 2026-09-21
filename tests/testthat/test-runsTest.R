

# -------------------------------------------------------------------------
# Basic functionality
# -------------------------------------------------------------------------
test_that("runsTest returns an htest object", {
  
  x <- c("S","S","T","S","T","T","T","S","T")
  
  res <- runsTest(x)
  
  expect_s3_class(res, "htest")
  expect_named(res$parameter, c("runs", "m", "n"))
  expect_true(is.numeric(res$p.value))
  expect_gte(res$p.value, 0)
  expect_lte(res$p.value, 1)
})
# -------------------------------------------------------------------------
# SPSS reference values
# -------------------------------------------------------------------------
test_that("SPSS example gives correct exact p-value", {
  
  x <- c(31,23,36,43,51,44,12,26,43,75,2,3,15,18,78,24,13,27,86,61,13,7,6,8)
  
  res <- runsTest(x, exact = TRUE)
  
  expect_equal(unname(res$parameter["runs"]), 10)
  expect_equal(unname(res$parameter["m"]),    12)
  expect_equal(unname(res$parameter["n"]),    12)
  expect_equal(res$p.value, 0.3009, tolerance = 1e-3)
})

test_that("SPSS example gives correct normal approximation", {
  
  x <- c(31,23,36,43,51,44,12,26,43,75,2,3,15,18,78,24,13,27,86,61,13,7,6,8)
  
  res <- runsTest(x, exact = FALSE)
  
  expect_equal(unname(res$statistic["z"]), -1.0436, tolerance = 1e-3)
  expect_equal(res$p.value, 0.2967, tolerance = 1e-3)
})

test_that("small example gives correct exact p-value", {
  
  x2 <- c(1,1,1,1,0,0,0,0,1,1)
  
  res <- runsTest(x2, exact = TRUE)
  
  expect_equal(res$p.value, 0.07143, tolerance = 1e-4)
})

test_that("small example gives correct normal approximation", {
  
  x2 <- c(1,1,1,1,0,0,0,0,1,1)
  
  res <- runsTest(x2, exact = FALSE)
  
  expect_equal(unname(res$statistic["z"]), -1.6156, tolerance = 1e-3)
  expect_equal(res$p.value, 0.1062, tolerance = 1e-3)
})
# -------------------------------------------------------------------------
# Median tie removal
# -------------------------------------------------------------------------
test_that("median removal gives same result as manual removal", {
  
  x <- c(13,3,14,14,1,14,3,8,14,17,9,14,13,2,16,1,3,12,13,14)
  
  res_auto <- runsTest(x)
  
  s      <- sign(x - median(x))
  x_bin  <- as.integer(s[s != 0] > 0)
  res_manual <- runsTest(x_bin)
  
  expect_equal(unname(res_auto$parameter), unname(res_manual$parameter))
  expect_equal(res_auto$p.value, res_manual$p.value)
})

test_that("all observations on one side after median removal throws error", {
  x <- c(1, 2, 3, 3, 3)   # median=3, alle nicht-median Werte unter Median
  expect_error(runsTest(x), "one side")
})

# -------------------------------------------------------------------------
# Category order preserved
# -------------------------------------------------------------------------
test_that("original category order is preserved for character input", {
  
  # "S" appears first -> coded as 0, "T" coded as 1
  x1 <- c("S","S","T","S","T","T","T","S","T")
  x2 <- c("T","T","S","T","S","S","S","T","S")   # reversed
  
  res1 <- runsTest(x1)
  res2 <- runsTest(x2)
  
  # runs and m/n swapped but p-value identical
  expect_equal(res1$p.value, res2$p.value)
  expect_equal(unname(res1$parameter["runs"]),
               unname(res2$parameter["runs"]))
})
# -------------------------------------------------------------------------
# Exact vs asymptotic auto-selection
# -------------------------------------------------------------------------
test_that("auto selects exact for n <= 30", {
  
  x <- c(1,1,1,1,0,0,0,0,1,1)
  
  res <- runsTest(x)
  
  expect_null(res$statistic)
})

test_that("auto selects normal for n > 30", {
  
  set.seed(1)
  x <- sample(0:1, 40, replace = TRUE)
  
  res <- runsTest(x)
  
  expect_false(is.null(res$statistic))
})
# -------------------------------------------------------------------------
# Wald-Wolfowitz two-sample test
# -------------------------------------------------------------------------
test_that("Wald-Wolfowitz gives correct exact p-value", {
  
  A <- c(35,44,39,50,48,29,60,75,49,66)
  B <- c(17,23,13,24,33,21,18,16,32)
  
  res <- runsTest(A, B, exact = TRUE)
  
  expect_equal(unname(res$parameter["runs"]), 4)
  expect_equal(res$p.value, 0.003139276, tolerance = 1e-5)
  expect_match(res$method, "Wald-Wolfowitz")
})

test_that("Wald-Wolfowitz gives correct normal approximation", {
  
  A <- c(35,44,39,50,48,29,60,75,49,66)
  B <- c(17,23,13,24,33,21,18,16,32)
  
  res <- runsTest(A, B, exact = FALSE)
  
  expect_equal(unname(res$statistic["z"]), -2.8287, tolerance = 1e-3)
  expect_equal(res$p.value, 0.004674, tolerance = 1e-4)
})

test_that("inter-group ties produce warning", {
  
  A <- c(10,14,17,19,34)
  B <- c(12,13,17,19,22)
  
  expect_warning(runsTest(A, B), "ties")
})

test_that("empty group throws error", {
  expect_error(runsTest(1:5, numeric(0)))
})
# -------------------------------------------------------------------------
# Alternatives
# -------------------------------------------------------------------------
test_that("less and greater p-values sum to 1 for normal", {
  
  x <- c(1,1,1,1,0,0,0,0,1,1)
  
  r_less    <- runsTest(x, alternative = "less",    exact = FALSE)
  r_greater <- runsTest(x, alternative = "greater", exact = FALSE)
  
  expect_equal(r_less$p.value + r_greater$p.value, 1, tolerance = 1e-10)
})

test_that("two-sided p >= one-sided p", {
  
  x <- c(1,1,1,1,0,0,0,0,1,1)
  
  r_less <- runsTest(x, alternative = "less",      exact = TRUE)
  r_two  <- runsTest(x, alternative = "two.sided", exact = TRUE)
  
  expect_gte(r_two$p.value, r_less$p.value)
})

test_that("invalid alternative throws error", {
  expect_error(runsTest(c(0,1,0,1), alternative = "invalid"))
})
# -------------------------------------------------------------------------
# NA handling
# -------------------------------------------------------------------------
test_that("NA with na.rm = TRUE gives same result as clean data", {
  
  x    <- c(1,0,1,1,0,0,1,0,1,1,0,1)
  x_na <- x
  x_na[c(3,8)] <- NA
  
  res_clean <- runsTest(x[!is.na(x_na)])
  res_na    <- runsTest(x_na, na.rm = TRUE)
  
  expect_equal(unname(res_clean$parameter), unname(res_na$parameter))
  expect_equal(res_clean$p.value, res_na$p.value)
})

test_that("NA with na.rm = FALSE throws error", {
  
  x <- c(1,0,NA,1,0)
  
  expect_error(runsTest(x, na.rm = FALSE), "NA")
})
# -------------------------------------------------------------------------
# Input validation
# -------------------------------------------------------------------------
test_that("single-value data throws error", {
  expect_error(runsTest(rep(1, 10)), "two distinct")
})

test_that("more than two categories throws error", {
  expect_error(runsTest(c("A","B","C","A","B")), "two distinct")
})

test_that("parameter values are correct", {
  
  x <- c(1,1,1,1,0,0,0,0,1,1)
  
  res <- runsTest(x, exact = TRUE)
  
  expect_equal(unname(res$parameter["runs"]), 3)
  expect_equal(unname(res$parameter["m"]),    4)
  expect_equal(unname(res$parameter["n"]),    6)
})
# -------------------------------------------------------------------------
# Print compatibility
# -------------------------------------------------------------------------
test_that("print.htest works", {
  
  x <- c(31,23,36,43,51,44,12,26,43,75,2,3,15,18,78,24,13,27,86,61,13,7,6,8)
  
  res <- runsTest(x)
  
  expect_output(print(res), "Runs")
})



# -- added --------------------------------------------------------------------

# exact null distribution of the number of runs (Swed & Eisenhart 1943)
druns <- function(n0, n1) {
  rr <- 2:(2 * min(n0, n1) + (n0 != n1))
  p <- vapply(rr, function(k) {
    j <- k %/% 2
    if (k %% 2 == 0) 2 * choose(n0 - 1, j - 1) * choose(n1 - 1, j - 1)
    else choose(n0 - 1, j) * choose(n1 - 1, j - 1) +
      choose(n0 - 1, j - 1) * choose(n1 - 1, j)
  }, 0) / choose(n0 + n1, n0)
  setNames(p, rr)
}

test_that("exact one-sided p-values are the tails of the runs distribution", {
  x <- c(1, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1)          # 7 runs, 6/7
  d <- druns(6, 7)
  r <- 7
  expect_equal(runsTest(x, exact = TRUE, alternative = "less")$p.value,
               sum(d[as.numeric(names(d)) <= r]))
  expect_equal(runsTest(x, exact = TRUE, alternative = "greater")$p.value,
               sum(d[as.numeric(names(d)) >= r]))
  expect_equal(sum(d), 1)
})

test_that("normal approximation by hand, with and without continuity correction", {
  x <- c(1, 1, 1, 1, 0, 0, 0, 0, 1, 1)
  n0 <- 4; n1 <- 6; R <- 3
  mu <- 1 + 2 * n0 * n1 / 10
  s2 <- 2 * n0 * n1 * (2 * n0 * n1 - 10) / (10^2 * 9)

  r <- runsTest(x, exact = FALSE, correct = FALSE)
  expect_equal(unname(r$statistic), (R - mu) / sqrt(s2))
  expect_equal(r$p.value, 2 * pnorm(-abs((R - mu) / sqrt(s2))))

  r <- runsTest(x, exact = FALSE)                          # R < mu: + 0.5
  expect_equal(unname(r$statistic), (R - mu + 0.5) / sqrt(s2))

  # R above mu: - 0.5
  y <- c(1, 0, 1, 0, 1, 0, 1, 0, 1, 1)
  mu <- 1 + 2 * 4 * 6 / 10
  expect_equal(unname(runsTest(y, exact = FALSE)$statistic),
               (9 - mu - 0.5) / sqrt(s2))
  # |R - mu| <= 0.5: no correction
  z <- c(1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1)               # R = 7, mu = 7
  expect_equal(unname(runsTest(z, exact = FALSE)$statistic),
               unname(runsTest(z, exact = FALSE, correct = FALSE)$statistic))
})

test_that("alternative texts", {
  x <- c(1, 1, 1, 1, 0, 0, 0, 0, 1, 1)
  expect_match(runsTest(x, alternative = "less")$alternative, "clustering")
  expect_match(runsTest(x, alternative = "greater")$alternative, "alternation")
  expect_match(runsTest(x)$alternative, "not equal")
})

test_that("logical and factor input", {
  x <- c(1, 1, 1, 1, 0, 0, 0, 0, 1, 1)
  ref <- runsTest(x)
  expect_equal(runsTest(x == 1)$p.value, ref$p.value)
  expect_equal(runsTest(factor(x, labels = c("a", "b")))$p.value, ref$p.value)
})

test_that("all NA and too few categories", {
  expect_error(runsTest(c(NA, NA)), "empty after removing NA")
  expect_error(runsTest(c("a", NA, "a")), "two distinct")
})

test_that("Wald-Wolfowitz equals the one-sample test on the merged group labels", {
  A <- c(35, 44, 39, 50, 48, 29, 60, 75, 49, 66)
  B <- c(17, 23, 13, 24, 33, 21, 18, 16, 32)
  lab <- c(rep(0, 10), rep(1, 9))[order(c(A, B))]
  ww <- runsTest(A, B, exact = FALSE, alternative = "less")
  one <- runsTest(lab, exact = FALSE, alternative = "less")
  expect_equal(ww$statistic, one$statistic)
  expect_equal(ww$p.value, one$p.value)
  expect_identical(ww$data.name, "A and B")
  expect_identical(unname(ww$parameter[c("m", "n")]), c(10L, 9L))
})

test_that("Wald-Wolfowitz removes missing values", {
  A <- c(35, 44, NA, 50, 48, 29, 60, 75, 49, 66)
  B <- c(17, 23, 13, NA, 33, 21, 18, 16, 32)
  r <- runsTest(A, B, exact = FALSE)
  # regression: NA was sorted last and counted as an observation
  expect_identical(unname(r$parameter[c("m", "n")]), c(9L, 8L))
  expect_equal(r$p.value, runsTest(A[!is.na(A)], B[!is.na(B)], exact = FALSE)$p.value)
  expect_error(runsTest(A, B, na.rm = FALSE), "NA")
  expect_error(runsTest(c(NA, NA), 1:3), "non-empty")
})

test_that("formula interface", {
  d <- data.frame(v = c(35, 44, 39, 50, 48, 29, 60, 75, 49, 66,
                        17, 23, 13, 24, 33, 21, 18, 16, 32),
                  g = factor(rep(c("A", "B"), c(10, 9))))
  f <- runsTest(v ~ g, data = d, exact = FALSE)
  ref <- runsTest(d$v[d$g == "A"], d$v[d$g == "B"], exact = FALSE)
  expect_equal(f$p.value, ref$p.value)
  expect_match(f$method, "Wald-Wolfowitz")
  # regression: data.name was read from a component that does not exist
  expect_type(f$data.name, "character")
  expect_match(f$data.name, "v")

  # regression: do.call() evaluated subset outside the data
  s <- runsTest(v ~ g, data = d, subset = v > 20, exact = FALSE)
  dd <- d[d$v > 20, ]
  expect_equal(s$p.value, runsTest(dd$v[dd$g == "A"], dd$v[dd$g == "B"],
                                   exact = FALSE)$p.value)

  expect_error(runsTest(~ v, data = d), "'formula' missing or incorrect")
})
