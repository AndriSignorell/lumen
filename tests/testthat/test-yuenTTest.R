
# -------------------------------------------------------------------------
# Basic functionality
# -------------------------------------------------------------------------

test_that("yuenTTest returns an htest object", {
  
  set.seed(1)
  
  x <- rnorm(20)
  y <- rnorm(20)
  
  res <- yuenTTest(x, y)
  
  expect_s3_class(res, "htest")
  
  expect_named(
    res$parameter,
    c("df", "trim")
  )
  
  expect_true(is.numeric(res$statistic))
  expect_true(is.numeric(res$p.value))
  
  expect_true(res$p.value >= 0)
  expect_true(res$p.value <= 1)
})


# -------------------------------------------------------------------------
# One-sample test
# -------------------------------------------------------------------------

test_that("one-sample test works", {
  
  set.seed(1)
  
  x <- rnorm(30, mean = 5)
  
  res <- yuenTTest(x, mu = 5)
  
  expect_s3_class(res, "htest")
  
  expect_equal(
    names(res$estimate),
    "trimmed mean of x"
  )
  
  expect_equal(
    res$null.value,
    c("trimmed mean" = 5)
  )
})


# -------------------------------------------------------------------------
# Two-sample test
# -------------------------------------------------------------------------

test_that("two-sample test works", {
  
  set.seed(1)
  
  x <- rnorm(25)
  y <- rnorm(30)
  
  res <- yuenTTest(x, y)
  
  expect_equal(
    names(res$estimate),
    c(
      "trimmed mean of x",
      "trimmed mean of y"
    )
  )
  
  expect_true(is.finite(res$statistic))
  expect_true(is.finite(res$parameter["df"]))
})


# -------------------------------------------------------------------------
# Paired test
# -------------------------------------------------------------------------

test_that("paired test uses paired differences", {
  
  set.seed(1)
  
  x <- rnorm(20)
  y <- x + rnorm(20, sd = 0.5)
  
  res <- yuenTTest(
    x,
    y,
    paired = TRUE
  )
  
  expect_match(
    res$method,
    "Paired"
  )
  
  expect_equal(
    names(res$estimate),
    "trimmed mean of the differences"
  )
  
  expect_true(is.finite(res$statistic))
})


# -------------------------------------------------------------------------
# Formula interface
# -------------------------------------------------------------------------

test_that("formula interface works", {
  
  set.seed(1)
  
  x <- rnorm(20)
  y <- rnorm(20)
  
  dat <- data.frame(
    value = c(x, y),
    group = factor(rep(
      c("A", "B"),
      c(length(x), length(y))
    ))
  )
  
  res1 <- yuenTTest(x, y)
  
  res2 <- yuenTTest(
    value ~ group,
    data = dat
  )
  
  expect_equal(
    unname(res1$statistic),
    unname(res2$statistic)
  )
  
  expect_equal(
    unname(res1$p.value),
    unname(res2$p.value)
  )
})


# -------------------------------------------------------------------------
# Alternative hypotheses
# -------------------------------------------------------------------------

test_that("alternative hypotheses work", {
  
  set.seed(1)
  
  x <- rnorm(20)
  y <- rnorm(20)
  
  res_less <- yuenTTest(
    x,
    y,
    alternative = "less"
  )
  
  res_greater <- yuenTTest(
    x,
    y,
    alternative = "greater"
  )
  
  res_two <- yuenTTest(
    x,
    y,
    alternative = "two.sided"
  )
  
  expect_true(
    all(is.finite(c(
      res_less$p.value,
      res_greater$p.value,
      res_two$p.value
    )))
  )
})


# -------------------------------------------------------------------------
# Trim validation
# -------------------------------------------------------------------------

test_that("invalid trim values throw errors", {
  
  x <- 1:10
  y <- 1:10
  
  expect_error(
    yuenTTest(x, y, trim = -0.1),
    "trim"
  )
  
  expect_error(
    yuenTTest(x, y, trim = 0.5),
    "trim"
  )
  
  expect_error(
    yuenTTest(x, y, trim = NA),
    "trim"
  )
  
  expect_error(
    yuenTTest(x, y, trim = c(0.1, 0.2)),
    "trim"
  )
})


# -------------------------------------------------------------------------
# Degrees of freedom protection
# -------------------------------------------------------------------------

test_that("too-large trim levels fail for small samples", {
  
  x <- 1:5
  y <- 1:5
  
  expect_error(
    yuenTTest(x, y, trim = 0.45),
    "trim level too large"
  )
})


# -------------------------------------------------------------------------
# Missing values
# -------------------------------------------------------------------------

test_that("missing values are removed", {
  
  x <- c(1, 2, 3, 4, NA, Inf)
  y <- c(1, 2, 3, 4, NaN)
  
  res <- yuenTTest(x, y)
  
  expect_true(is.finite(res$statistic))
  expect_true(is.finite(res$p.value))
})


# -------------------------------------------------------------------------
# Constant data
# -------------------------------------------------------------------------

test_that("constant data throw errors", {
  
  x <- rep(1, 20)
  y <- rep(1, 20)
  
  expect_error(
    yuenTTest(x),
    "essentially constant"
  )
  
  expect_error(
    yuenTTest(x, y),
    "essentially constant"
  )
})


# -------------------------------------------------------------------------
# Confidence intervals
# -------------------------------------------------------------------------

test_that("confidence intervals are returned correctly", {
  
  set.seed(1)
  
  x <- rnorm(20)
  y <- rnorm(20)
  
  res <- yuenTTest(
    x,
    y,
    conf.level = 0.90
  )
  
  expect_equal(
    length(res$conf.int),
    2
  )
  
  expect_equal(
    attr(res$conf.int, "conf.level"),
    0.90
  )
  
  expect_named(
    res$conf.int,
    c("lower", "upper")
  )
})


# -------------------------------------------------------------------------
# Parameter validation
# -------------------------------------------------------------------------

test_that("invalid mu values throw errors", {
  
  x <- rnorm(20)
  
  expect_error(
    yuenTTest(x, mu = c(1, 2)),
    "mu"
  )
  
  expect_error(
    yuenTTest(x, mu = NA),
    "mu"
  )
})


test_that("invalid conf.level values throw errors", {
  
  x <- rnorm(20)
  
  expect_error(
    yuenTTest(x, conf.level = 2),
    "conf.level"
  )
  
  expect_error(
    yuenTTest(x, conf.level = 0),
    "conf.level"
  )
})


# -------------------------------------------------------------------------
# Paired test validation
# -------------------------------------------------------------------------

test_that("paired test requires y", {
  
  x <- rnorm(20)
  
  expect_error(
    yuenTTest(
      x,
      paired = TRUE
    ),
    "'y' is missing"
  )
})


# -------------------------------------------------------------------------
# Numerical sanity
# -------------------------------------------------------------------------

test_that("p-values remain in [0,1]", {
  
  set.seed(1)
  
  for(i in 1:100) {
    
    x <- rnorm(sample(10:50, 1))
    y <- rnorm(sample(10:50, 1))
    
    res <- yuenTTest(x, y)
    
    expect_true(is.finite(res$p.value))
    
    expect_true(res$p.value >= 0)
    expect_true(res$p.value <= 1)
  }
})


# -------------------------------------------------------------------------
# Print compatibility
# -------------------------------------------------------------------------

test_that("print.htest works", {
  
  set.seed(1)
  
  x <- rnorm(20)
  y <- rnorm(20)
  
  res <- yuenTTest(x, y)
  
  expect_output(
    print(res),
    "Yuen"
  )
})


test_that("formula interface splits by group, not x vs convenience-y", {
  
  set.seed(1)
  
  x <- rnorm(20)
  y <- rnorm(20, mean = 3)   # deutlich verschieden, damit ein Bug auffällt
  
  dat <- data.frame(
    value = c(x, y),
    group = factor(rep(c("A", "B"), c(20, 20)))
  )
  
  res1 <- yuenTTest(x, y)
  res2 <- yuenTTest(value ~ group, data = dat)
  
  # regression guard: resolveFormula()'s 'x' is now the full response
  # (both groups), not group A alone - yuenTTest.formula() must split
  # by 'group' explicitly, or this compares "all 40 obs" against
  # "group B only" and silently returns the wrong statistic.
  expect_equal(unname(res1$statistic), unname(res2$statistic))
  expect_equal(unname(res1$p.value), unname(res2$p.value))
})



test_that("two-sample results agree with WRS / WRS2 / PairedData (issue #47)", {
  
  res <- yuenTTest(extra ~ group, data = sleep)
  
  # reference: WRS::yuen(), WRS2::yuen(), PairedData::yuen.t.test()
  expect_equal(unname(res$statistic), -1.616777, tolerance = 1e-6)
  expect_equal(unname(res$parameter["df"]), 8.264709, tolerance = 1e-6)
  expect_equal(res$p.value, 0.1433783, tolerance = 1e-6)
  expect_equal(unname(res$conf.int[1:2]), c(-4.0306400, 0.6973066),
               tolerance = 1e-6)
})


test_that("the confidence interval is for the parameter, not shifted by mu", {
  
  x <- sleep$extra[sleep$group == 1]
  y <- sleep$extra[sleep$group == 2]
  
  # two-sample: centred at the difference of the trimmed means
  ci <- yuenTTest(x, y)$conf.int
  expect_equal(mean(ci), mean(x, trim = 0.2) - mean(y, trim = 0.2))
  
  # independent of mu, for all three designs
  expect_equal(yuenTTest(x, y, mu = 1)$conf.int, ci)
  expect_equal(yuenTTest(x, mu = 1)$conf.int, yuenTTest(x)$conf.int)
  expect_equal(yuenTTest(x, y, paired = TRUE, mu = 1)$conf.int,
               yuenTTest(x, y, paired = TRUE)$conf.int)
  
  # one-sided bounds are consistent with the two-sided interval
  lo <- yuenTTest(x, y, alternative = "greater", conf.level = 0.975)$conf.int
  up <- yuenTTest(x, y, alternative = "less",    conf.level = 0.975)$conf.int
  expect_equal(unname(c(lo[1], up[2])), unname(ci[1:2]))
})


test_that("winsorizing is done at order statistics", {
  
  # n = 10, trim = 0.2: g = 2, so the two smallest become z[3] and the two
  # largest z[8], whatever quantile() would interpolate
  z <- c(-10, -5, 1, 2, 3, 4, 5, 6, 20, 30)
  w <- c(1, 1, 1, 2, 3, 4, 5, 6, 6, 6)
  
  expect_equal(.winsorVar(z, 0.2), var(w))
  # n = 10, g = 2, h = 6
  expect_equal(.trimmedSE(z, 0.2), sqrt(9 * var(w) / (6 * 5)))
  
  # no trimming: plain variance
  expect_equal(.winsorVar(z, 0), var(z))
})


test_that("constant data are rejected, also when all values are 0", {
  
  msg <- "essentially constant"
  expect_error(yuenTTest(rep(0, 20)), msg)
  expect_error(yuenTTest(rep(0, 20), rep(0, 20)), msg)
  expect_error(yuenTTest(1:20, 1:20, paired = TRUE), msg)
})


test_that("without trimming the tests are the ordinary t-tests", {
  
  set.seed(47)
  x <- rnorm(15)
  y <- rnorm(12, 1, 2)
  xp <- x[1:12]
  
  same <- function(a, b) {
    expect_equal(unname(a$statistic), unname(b$statistic))
    expect_equal(unname(a$parameter["df"]), unname(b$parameter["df"]))
    expect_equal(a$p.value, b$p.value)
    expect_equal(unname(a$conf.int[1:2]), unname(b$conf.int[1:2]))
  }
  
  same(yuenTTest(x, mu = 0.3, trim = 0), t.test(x, mu = 0.3))
  same(yuenTTest(x, y, trim = 0),         t.test(x, y))
  same(yuenTTest(xp, y, paired = TRUE, trim = 0),
       t.test(xp, y, paired = TRUE))
  
  # nothing trimmed although trim > 0: n = 4, trim = 0.2 gives g = 0
  same(yuenTTest(c(1, 4, 2, 7), trim = 0.2), t.test(c(1, 4, 2, 7)))
})


test_that("paired equals the one-sample test on complete differences", {
  
  x <- c(5.1, NA, 4.8, 6.0, 5.5, 7.2, 4.9, 5.8, 6.3, 5.0, 4.4, 6.8)
  y <- c(4.9, 5.0, NA, 5.1, 5.6, 6.0, 4.1, 5.9, 5.2, 4.2, 4.5, 5.9)
  
  d <- (x - y)[complete.cases(x, y)]
  a <- yuenTTest(x, y, paired = TRUE, mu = 0.1)
  b <- yuenTTest(d, mu = 0.1)
  
  expect_equal(a$statistic, b$statistic)
  expect_equal(a$parameter, b$parameter)
  expect_equal(a$p.value, b$p.value)
  expect_equal(a$conf.int, b$conf.int)
})


test_that("argument checks", {
  
  expect_error(yuenTTest(extra ~ group, data = sleep, paired = TRUE),
               "formula interface")
  # abbreviated: would be matched partially by yuenTTest.default()
  expect_error(yuenTTest(extra ~ group, data = sleep, pair = TRUE),
               "formula interface")
  expect_error(yuenTTest(extra ~ group, data = sleep, paired = NA),
               "formula interface")
  expect_s3_class(yuenTTest(extra ~ group, data = sleep, paired = FALSE),
                  "htest")

  expect_error(yuenTTest(1:10, 2:11, paired = NA), "'paired'")
  expect_error(yuenTTest(1:10, 2:11, paired = c(TRUE, FALSE)), "'paired'")
  expect_error(yuenTTest(1:10, mu = Inf), "'mu'")
})


test_that("formula interface honours subset", {
  # yuenTTest.formula() passes the subset through do.call(); zTest() had a
  # regression there, where the expression was evaluated outside the data
  f <- yuenTTest(extra ~ group, data = sleep, subset = ID != "1")
  s <- sleep[sleep$ID != "1", ]
  d <- yuenTTest(s$extra[s$group == 1], s$extra[s$group == 2])
  expect_equal(f$statistic, d$statistic)
  expect_equal(f$conf.int, d$conf.int)
})
