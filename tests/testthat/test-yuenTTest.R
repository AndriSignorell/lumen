
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
    c("trimmed mean difference" = 5)
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
    "difference in trimmed means"
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


