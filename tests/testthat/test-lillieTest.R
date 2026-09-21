library(testthat)
library(lumen)

test_that("lillieTest: returns htest", {
  expect_s3_class(lillieTest(rnorm(50)), "htest")
})

test_that("lillieTest: p.value in [0,1]", {
  set.seed(1)
  res <- lillieTest(rnorm(50))
  expect_true(res$p.value >= 0 && res$p.value <= 1)
})

test_that("lillieTest: normal data gives large p", {
  set.seed(42)
  expect_gt(lillieTest(rnorm(200))$p.value, 0.05)
})

test_that("lillieTest: non-normal data gives small p", {
  set.seed(1)
  expect_lt(lillieTest(runif(200))$p.value, 0.05)
})

test_that("lillieTest: statistic named K", {
  res <- lillieTest(rnorm(20))
  expect_named(res$statistic, "D")
})

test_that("lillieTest: statistic > 0", {
  set.seed(1)
  expect_gt(lillieTest(rnorm(30))$statistic, 0)
})

test_that("lillieTest: n < 5 throws error", {
  expect_error(lillieTest(rnorm(4)))
})

test_that("lillieTest: NAs silently removed", {
  set.seed(1)
  x <- c(rnorm(50), NA, NA)
  expect_s3_class(lillieTest(x), "htest")
})


# -- added --------------------------------------------------------------------

test_that("lillieTest: identical to nortest::lillie.test in every p-value branch", {
  skip_if_not_installed("nortest")
  set.seed(11)
  cases <- list(
    small_normal = rnorm(20),            # p > 0.1, KK <= 0.302 or low branch
    mid          = rt(40, 5),
    skewed       = rexp(60),             # small p, exp formula
    large_n      = rnorm(300),           # n > 100: Kd, nd = 100
    large_skew   = rexp(250),
    uniform      = runif(80),
    tiny         = c(1.2, 3.4, 2.2, 5.9, 4.1)
  )
  # a few extra draws to reach the polynomial branches of the p > 0.1 part
  for (i in 1:40) cases[[paste0("r", i)]] <- rnorm(25 + i)
  for (nm in names(cases)) {
    a <- lillieTest(cases[[nm]])
    b <- nortest::lillie.test(cases[[nm]])
    expect_equal(unname(a$statistic), unname(b$statistic), info = nm)
    expect_equal(a$p.value, b$p.value, info = nm)
  }
})

test_that("lillieTest: statistic is the KS distance to the fitted normal", {
  set.seed(2)
  x <- rexp(30)
  z <- (x - mean(x)) / sd(x)
  expect_equal(unname(lillieTest(x)$statistic),
               unname(suppressWarnings(ks.test(z, "pnorm"))$statistic))
})

test_that("lillieTest: invariant under location and scale", {
  set.seed(3)
  x <- rgamma(40, 2)
  expect_equal(lillieTest(x)$p.value, lillieTest(10 + 3 * x)$p.value)
})

test_that("lillieTest: input checks", {
  expect_error(lillieTest(letters), "numeric")
  expect_error(lillieTest(rep(2, 10)), "identical")
  expect_error(lillieTest(c(1:10, Inf)), "infinite")
  expect_error(lillieTest(c(1:4, NA, NA)), "at least 5")
  expect_identical(lillieTest(c(1, 4, 2, 8, 5, 7))$data.name,
                   "c(1, 4, 2, 8, 5, 7)")
})
