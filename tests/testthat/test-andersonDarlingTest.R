# Anderson-Darling test: statistic, htest contract, Braun adjustment,
# helper functions. Regression against goftest lives in test-ad-regression.R.

# textbook statistic, independent of the C++ code
adStatRef <- function(u) {
  u <- sort(u)
  n <- length(u)
  i <- seq_len(n)
  -n - sum((2 * i - 1) * (log(u) + log1p(-rev(u)))) / n
}


test_that("statistic matches the textbook formula", {
  set.seed(11)
  for (n in c(1, 2, 7, 50, 400)) {
    u <- runif(n)
    expect_equal(unname(andersonDarlingTest(u)$statistic), adStatRef(u),
                 tolerance = 1e-12)
  }
  u <- runif(25)
  expect_equal(ad_stat_cpp(sort(u)), adStatRef(u), tolerance = 1e-12)
  r <- ad_test_r_cpp(sort(u))
  expect_named(r, c("adstat", "pvalue"))
  expect_equal(r$pvalue, 1 - ad_prob_n_cpp(r$adstat, 25L), tolerance = 1e-15)
})


test_that("htest contract and method text", {
  set.seed(1)
  x <- rexp(20, rate = 0.5)
  r <- andersonDarlingTest(x, "pexp", rate = 0.5)

  expect_s3_class(r, "htest")
  expect_named(r, c("statistic", "p.value", "method", "data.name"))
  expect_identical(names(r$statistic), "An")
  expect_identical(r$data.name, "x")
  expect_true(r$p.value >= 0 && r$p.value <= 1)
  expect_match(r$method, "null hypothesis: exponential distribution", fixed = TRUE)
  expect_match(r$method, "with parameter rate = 0.5", fixed = TRUE)
  expect_match(r$method, "parameters fixed", fixed = TRUE)

  r2 <- andersonDarlingTest(x, "pgamma", shape = 1, rate = 0.5)
  expect_match(r2$method, "with parameters shape = 1, rate = 0.5", fixed = TRUE)
  # same null distribution, same statistic
  expect_equal(r2$statistic, r$statistic, tolerance = 1e-12)
})


test_that("null: name, name without leading 'p', function object, anonymous function", {
  set.seed(2)
  x <- rnorm(15)
  ref <- andersonDarlingTest(x, "pnorm")

  expect_equal(andersonDarlingTest(x, "norm")$statistic, ref$statistic)
  expect_equal(andersonDarlingTest(x, pnorm)$statistic, ref$statistic)
  expect_match(andersonDarlingTest(x, pnorm)$method, "Normal distribution",
               fixed = TRUE)

  r <- andersonDarlingTest(x, function(q) pnorm(q))
  expect_equal(r$statistic, ref$statistic)
  expect_match(r$method, "distribution .function\\(q\\) pnorm\\(q\\).")

  expect_match(andersonDarlingTest(x, "pnorm", nullname = "my H0")$method,
               "null hypothesis: my H0;", fixed = TRUE)
})


test_that("data.name and null name stay single strings for long expressions", {
  r <- andersonDarlingTest(
    c(0.05, 0.12, 0.18, 0.27, 0.33, 0.41, 0.46, 0.52, 0.58, 0.66, 0.71,
      0.79, 0.84, 0.91, 0.97),
    function(q) {
      z <- q
      punif(z)
    })
  expect_length(r$data.name, 1L)
  expect_length(gregexpr("null hypothesis:", r$method, fixed = TRUE)[[1L]], 1L)
})


test_that("missing values are removed", {
  set.seed(3)
  u <- runif(30)
  r1 <- andersonDarlingTest(u)
  r2 <- andersonDarlingTest(c(u[1:10], NA, u[11:30], NA))
  expect_equal(r2$statistic, r1$statistic)
  expect_equal(r2$p.value, r1$p.value)
})


test_that("p-value stays in [0, 1] for near-perfect fits", {
  # Marsaglia's finite-n correction goes below 0 for tiny statistics
  # (n = 5: A2 < 0.145), which used to give p = 1.000265 here
  for (n in 2:30) {
    p <- andersonDarlingTest((seq_len(n) - 0.5) / n)$p.value
    expect_true(p >= 0 && p <= 1, label = paste("p-value for n =", n))
  }
})


test_that("observations on the boundary of the support give An = Inf, p = 0", {
  r <- andersonDarlingTest(c(0, 0.3, 0.6))
  expect_identical(unname(r$statistic), Inf)
  expect_identical(r$p.value, 0)
})


test_that("estimated = TRUE uses Braun's adjustment", {
  set.seed(3)
  x <- rnorm(30)
  n <- length(x)
  m <- round(sqrt(n))

  set.seed(99)
  r <- andersonDarlingTest(x, "pnorm", mean = mean(x), sd = sd(x),
                           estimated = TRUE)

  expect_identical(names(r$statistic), "Anmax")
  expect_match(r$method, paste("Braun's adjustment using", m, "groups"),
               fixed = TRUE)
  expect_match(r$method, "parameters estimated from data", fixed = TRUE)

  # reconstruct with the same random grouping
  set.seed(99)
  g  <- factor(sample(seq_len(n) %% m))
  uu <- split(pnorm(x, mean(x), sd(x)), g)
  pg <- vapply(uu, function(u) .simpleADtest(u)$pvalue, 0)

  expect_equal(unname(r$statistic), max(vapply(uu, adStatRef, 0)),
               tolerance = 1e-12)
  expect_equal(r$p.value, 1 - (1 - min(pg))^m, tolerance = 1e-12)
})


test_that("estimated = TRUE with n <= 4 falls back to the simple test", {
  x <- c(-1, 0.2, 0.5, 1.3)
  expect_warning(r <- andersonDarlingTest(x, "pnorm", estimated = TRUE),
                 "too few observations")
  expect_identical(names(r$statistic), "An")
  expect_match(r$method, "parameters fixed", fixed = TRUE)
  expect_equal(r$statistic, andersonDarlingTest(x, "pnorm")$statistic)
})


test_that(".braun() refuses too few observations per group", {
  expect_error(.braun(runif(5), .simpleADtest, m = 3), "Insufficient data")
})


test_that("input errors", {
  expect_error(andersonDarlingTest(letters), "numeric")
  expect_error(andersonDarlingTest(c(NA_real_, NA_real_)), "not enough")
  expect_error(andersonDarlingTest(runif(5), function(q) 2 * q),
               "outside \\[0,1\\]")
  expect_error(andersonDarlingTest(runif(5), "doesNotExist123"),
               "should be a function")
})


test_that(".recogniseCdf()", {
  expect_identical(.recogniseCdf("pnorm"), "Normal distribution")
  expect_identical(.recogniseCdf("t"), "Student's t distribution")
  expect_identical(.recogniseCdf("unif"), "uniform distribution")
  expect_identical(.recogniseCdf("pAD"),
                   "null distribution of Anderson-Darling Test Statistic")
  expect_identical(.recogniseCdf("pCvM"),
                   "null distribution of Cramer-von Mises Test Statistic")
  expect_null(.recogniseCdf("pfoo"))
  expect_null(.recogniseCdf("p"))
  expect_null(.recogniseCdf(""))
  expect_null(.recogniseCdf(1))
  expect_null(.recogniseCdf(c("pnorm", "pexp")))
})


test_that(".getCdf()", {
  expect_identical(.getCdf("t"), stats::pt)          # leading 'p' added
  expect_identical(.getCdf("pexp"), stats::pexp)
  expect_identical(.getCdf(pnorm), pnorm)
  expect_null(.getCdf("doesNotExist123", fatal = FALSE))
  expect_null(.getCdf(42, fatal = FALSE))
  expect_error(.getCdf(42), "should be a function")
  expect_null(.getfunky("doesNotExist123"))
})
