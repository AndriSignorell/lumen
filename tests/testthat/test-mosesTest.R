
# -------------------------------------------------------------------------
# Basic functionality
# -------------------------------------------------------------------------

test_that("mosesTest returns an htest-compatible object", {
  
  x <- c(0.80, 0.83, 1.89, 1.04, 1.45,
         1.38, 1.91, 1.64, 0.73, 1.46)
  
  y <- c(1.15, 0.88, 0.90, 0.74, 1.21)
  
  res <- mosesTest(x, y)
  
  expect_s3_class(res, "mosesTestResult")
  expect_s3_class(res, "htest")
  
  expect_named(
    res$statistic,
    c("S_raw", "S_trimmed")
  )
  
  expect_named(
    res$p.value,
    c("p_raw", "p_trimmed")
  )
  
  expect_true(all(res$p.value >= 0))
  expect_true(all(res$p.value <= 1))
})


# -------------------------------------------------------------------------
# Formula interface
# -------------------------------------------------------------------------

test_that("formula interface works correctly", {
  
  x <- c(1, 2, 3, 10, 20)
  y <- c(2, 2, 3, 4)
  
  df <- data.frame(
    value = c(x, y),
    group = factor(rep(
      c("control", "experiment"),
      c(length(x), length(y))
    ))
  )
  
  res1 <- mosesTest(x, y)
  
  res2 <- mosesTest(
    value ~ group,
    data = df
  )
  
  expect_equal(
    res1$statistic,
    res2$statistic
  )
  
  expect_equal(
    res1$p.value,
    res2$p.value
  )
})


# -------------------------------------------------------------------------
# Extreme parameter validation
# -------------------------------------------------------------------------

test_that("invalid extreme values throw errors", {
  
  x <- 1:10
  y <- 1:5
  
  expect_error(
    mosesTest(x, y, extreme = -1),
    "non-negative integer"
  )
  
  expect_error(
    mosesTest(x, y, extreme = 1.5),
    "non-negative integer"
  )
  
  expect_error(
    mosesTest(x, y, extreme = NA),
    "non-negative integer"
  )
  
  expect_error(
    mosesTest(x, y, extreme = c(1, 2)),
    "non-negative integer"
  )
})


# -------------------------------------------------------------------------
# Excessive trimming reduced
# -------------------------------------------------------------------------

test_that("extreme values larger than allowed are reduced", {
  
  x <- 1:6
  y <- 1:4
  
  expect_warning(
    
    res <- mosesTest(
      x,
      y,
      extreme = 100
    ),
    
    "reduced"
  )
  
  expect_equal(
    res$extreme,
    floor((length(x) - 2L) / 2L)
  )
})


# -------------------------------------------------------------------------
# Tie handling
# -------------------------------------------------------------------------

test_that("ties.method argument is validated", {
  
  x <- c(1, 1, 2, 3, 4)
  y <- c(1, 2, 2)
  
  expect_error(
    mosesTest(
      x,
      y,
      ties.method = "invalid"
    )
  )
})


test_that("approximate p-value warning appears for non-first ties", {
  
  x <- c(1, 1, 2, 3, 4)
  y <- c(1, 2, 2)
  
  expect_warning(
    
    mosesTest(
      x,
      y,
      ties.method = "average"
    ),
    
    "approximate p-values"
  )
})


test_that("ties.method = 'first' is reproducible", {
  
  x <- c(1, 1, 2, 3, 4)
  y <- c(1, 2, 2)
  
  res1 <- mosesTest(
    x,
    y,
    ties.method = "first"
  )
  
  res2 <- mosesTest(
    x,
    y,
    ties.method = "first"
  )
  
  expect_equal(
    res1$statistic,
    res2$statistic
  )
  
  expect_equal(
    res1$p.value,
    res2$p.value
  )
})


# -------------------------------------------------------------------------
# Numerical properties
# -------------------------------------------------------------------------

test_that("p-values remain within [0,1]", {
  
  set.seed(1)
  
  for (i in 1:100) {
    
    x <- rnorm(sample(5:30, 1))
    y <- rnorm(sample(5:30, 1))
    
    res <- mosesTest(x, y)
    
    expect_true(all(is.finite(res$p.value)))
    
    expect_true(all(res$p.value >= 0))
    expect_true(all(res$p.value <= 1))
  }
})


# -------------------------------------------------------------------------
# Trimming changes statistic
# -------------------------------------------------------------------------

test_that("trimming affects span statistic", {
  
  x <- c(1, 2, 3, 4, 100, 101, 102)
  y <- c(5, 6, 7, 8)
  
  res0 <- mosesTest(
    x,
    y,
    extreme = 0
  )
  
  res1 <- mosesTest(
    x,
    y,
    extreme = 1
  )
  
  expect_true(
    res1$statistic["S_trimmed"] <=
      res0$statistic["S_raw"]
  )
})


# -------------------------------------------------------------------------
# Missing and infinite values
# -------------------------------------------------------------------------

test_that("non-finite observations are removed", {
  
  x <- c(1, 2, 3, 4, NA, Inf)
  y <- c(1, 2, 3, NaN)
  
  res <- mosesTest(x, y)
  
  expect_equal(
    unname(res$parameter["n_control"]),
    4
  )
  
  expect_equal(
    unname(res$parameter["n_experiment"]),
    3
  )
})


# -------------------------------------------------------------------------
# Small sample protection
# -------------------------------------------------------------------------

test_that("too-small samples throw errors", {
  
  expect_error(
    mosesTest(1:3, 1:5),
    "at least 4"
  )
  
  expect_error(
    mosesTest(1:5, numeric(0)),
    "at least 1"
  )
})


# -------------------------------------------------------------------------
# Print method
# -------------------------------------------------------------------------

test_that("print method works", {
  
  x <- 1:10
  y <- 1:5
  
  res <- mosesTest(x, y)
  
  expect_output(
    print(res),
    "Moses Test of Extreme Reactions"
  )
  
  expect_output(
    print(res),
    "Without trimming"
  )
})