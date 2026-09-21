library(testthat)
library(lumen)

# Equal variances: homoscedastic groups
set.seed(1)
g1 <- rnorm(30, 0, 1); g2 <- rnorm(30, 5, 1); g3 <- rnorm(30, 10, 1)
df_eq <- data.frame(x = c(g1, g2, g3), g = factor(rep(1:3, each = 30)))

# Unequal variances: heteroscedastic groups
g4 <- rnorm(30, 0, 1); g5 <- rnorm(30, 0, 5)
df_uneq <- data.frame(x = c(g4, g5), g = factor(rep(1:2, each = 30)))

test_that("leveneTest: returns htest (formula method)", {
  expect_s3_class(leveneTest(x ~ g, data = df_eq), "htest")
})

test_that("leveneTest: p.value in [0,1]", {
  res <- leveneTest(x ~ g, data = df_eq)
  expect_true(res$p.value >= 0 && res$p.value <= 1)
})

test_that("leveneTest: equal variances gives large p", {
  set.seed(42)
  g1 <- rnorm(50, sd=1); g2 <- rnorm(50, sd=1)
  df <- data.frame(x = c(g1,g2), g = factor(rep(1:2, each=50)))
  expect_gt(leveneTest(x ~ g, data = df)$p.value, 0.05)
})

test_that("leveneTest: unequal variances gives small p", {
  set.seed(42)
  g1 <- rnorm(100, sd=1); g2 <- rnorm(100, sd=10)
  df <- data.frame(x = c(g1,g2), g = factor(rep(1:2, each=100)))
  expect_lt(leveneTest(x ~ g, data = df)$p.value, 0.05)
})

test_that("leveneTest: center=mean gives original Levene test", {
  res <- leveneTest(x ~ g, data = df_eq, center = mean)
  expect_s3_class(res, "htest")
})

test_that("leveneTest: default method (formula) works", {
  res <- leveneTest(x ~ g, data = df_uneq)
  expect_false(is.null(res$statistic))
})

test_that("leveneTest: default method (default) works with vector + factor", {
  res <- leveneTest(df_eq$x, df_eq$g)
  expect_s3_class(res, "htest")
})


# leveneTest() ---------------------------------------------------------------

d <- InsectSprays

manualF <- function(x, g, center = median, ...) {
  g <- factor(g)
  dev <- abs(x - ave(x, g, FUN = function(z) center(z, ...)))
  a <- anova(lm(dev ~ g))
  list(F = a$`F value`[1], df = a$Df)
}

test_that("statistic and df agree with the ANOVA of absolute deviations", {
  r <- leveneTest(d$count, d$spray)
  m <- manualF(d$count, d$spray)
  
  expect_s3_class(r, "htest")
  expect_equal(r$statistic, c(F = m$F))
  expect_equal(r$parameter, c("num df" = 5, "denom df" = 66))
  expect_equal(unname(r$parameter), m$df)
  expect_equal(r$p.value, pf(m$F, 5, 66, lower.tail = FALSE))
  expect_s3_class(r$anova_tab, "anova")
  expect_match(r$method, "center = median", fixed = TRUE)
})

test_that("regression value: count ~ spray, InsectSprays", {
  expect_equal(unname(leveneTest(count ~ spray, data = d)$statistic),
               3.821356, tolerance = 1e-6)
})

test_that("formula, default and list input agree", {
  f <- leveneTest(count ~ spray, data = d)
  g <- leveneTest(d$count, d$spray)
  l <- leveneTest(split(d$count, d$spray))
  
  expect_equal(f$statistic, g$statistic)
  expect_equal(l$statistic, g$statistic)
  expect_equal(l$parameter, g$parameter)
  expect_match(f$data.name, "count")
})

test_that("center = mean and extra arguments to center", {
  r <- leveneTest(count ~ spray, data = d, center = mean)
  expect_equal(unname(r$statistic), manualF(d$count, d$spray, mean)$F)
  expect_match(r$method, "center = mean)", fixed = TRUE)
  
  r <- leveneTest(count ~ spray, data = d, center = mean, trim = 0.1)
  expect_equal(unname(r$statistic),
               manualF(d$count, d$spray, mean, trim = 0.1)$F)
  expect_match(r$method, "center = mean(trim=0.1)", fixed = TRUE)
  
  r <- leveneTest(d$count, d$spray, center = mean, trim = 0.2)
  expect_match(r$method, "mean(trim=0.2)", fixed = TRUE)
})

test_that("formula subset", {
  r <- leveneTest(count ~ spray, data = d, subset = spray %in% c("A", "B", "C"))
  s <- d[d$spray %in% c("A", "B", "C"), ]
  expect_equal(r$statistic, leveneTest(s$count, droplevels(s$spray))$statistic)
  expect_equal(unname(r$parameter[1]), 2)
})

test_that("missing values are dropped", {
  x <- c(d$count, NA, 5)
  g <- factor(c(as.character(d$spray), "A", NA))
  expect_equal(leveneTest(x, g)$statistic, leveneTest(d$count, d$spray)$statistic)
  
  l <- split(d$count, d$spray); l$A <- c(l$A, NA)
  expect_equal(leveneTest(l)$statistic, leveneTest(d$count, d$spray)$statistic)
})

test_that("list input: checks and warnings", {
  l <- split(d$count, d$spray)
  expect_warning(leveneTest(l, g = d$spray), "ignoring argument 'g'")
  expect_error(leveneTest(l[1]), "at least 2 elements")
  expect_error(leveneTest(list(1:3, numeric(0))), "all groups must contain data")
  expect_error(leveneTest(list(1:3, c(NA_real_, NA_real_))), "all groups must contain data")
  expect_warning(leveneTest(list(c(1, 4, 2), c(TRUE, FALSE, TRUE))),
                 "not numeric")
})

test_that("default method: checks", {
  expect_error(leveneTest(1:5, factor(c(1, 1, 2, 2))), "same length")
  expect_error(leveneTest(1:4, factor(rep("a", 4))), "same group")
  expect_error(leveneTest(1:4, factor(c("a", "a", NA, NA))), "same group")
})
