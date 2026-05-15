
# -------------------------------------------------------------------------
# Basic functionality
# -------------------------------------------------------------------------
test_that("vonNeumannTest returns an htest object", {
  
  set.seed(2)
  x <- runif(20)
  
  res <- vonNeumannTest(x)
  
  expect_s3_class(res, "htest")
  expect_named(res$statistic, "z")
  expect_named(res$parameter, "n")
  expect_true(is.numeric(res$p.value))
  expect_gte(res$p.value, 0)
  expect_lte(res$p.value, 1)
  expect_true(!is.null(res$vn))
})
# -------------------------------------------------------------------------
# Known reference values
# -------------------------------------------------------------------------
test_that("runif example gives correct z and p-value", {
  
  set.seed(2)
  x <- runif(20)
  
  res <- vonNeumannTest(x)
  
  expect_equal(unname(res$statistic["z"]), 0.73709, tolerance = 1e-4)
  expect_equal(res$p.value, 0.4611, tolerance = 1e-3)
})

test_that("tourists example gives correct z and p-value (less)", {
  
  tourists <- c(12362,12739,13057,13955,14123,15698,17523,18610,19842,
                20310,22500,23080,21916)
  
  res <- vonNeumannTest(tourists, alternative = "less")
  
  expect_equal(unname(res$statistic["z"]), -3.7369, tolerance = 1e-3)
  expect_equal(res$p.value, 0.00009316035, tolerance = 1e-6)
})
# -------------------------------------------------------------------------
# unbiased vs biased
# -------------------------------------------------------------------------
test_that("unbiased and biased give same z-statistic", {
  
  set.seed(1)
  x <- rnorm(30)
  
  r_unbiased <- vonNeumannTest(x, unbiased = TRUE)
  r_biased   <- vonNeumannTest(x, unbiased = FALSE)
  
  expect_equal(
    unname(r_unbiased$statistic["z"]),
    unname(r_biased$statistic["z"]),
    tolerance = 1e-10
  )
})

test_that("unbiased VN differs from biased VN", {
  
  set.seed(1)
  x <- rnorm(30)
  
  r_unbiased <- vonNeumannTest(x, unbiased = TRUE)
  r_biased   <- vonNeumannTest(x, unbiased = FALSE)
  
  expect_false(
    isTRUE(all.equal(r_unbiased$vn, r_biased$vn))
  )
})
# -------------------------------------------------------------------------
# Alternatives
# -------------------------------------------------------------------------
test_that("less alternative detects trend", {
  
  x <- cumsum(rnorm(30, mean = 1))
  
  res <- vonNeumannTest(x, alternative = "less")
  
  expect_lt(res$p.value, 0.05)
})

test_that("greater alternative detects oscillation", {
  
  x <- rep(c(-5, 5), 15) + rnorm(30, sd = 0.1)
  
  res <- vonNeumannTest(x, alternative = "greater")
  
  expect_lt(res$p.value, 0.05)
})

test_that("two-sided p-value >= one-sided p-value", {
  
  tourists <- c(12362,12739,13057,13955,14123,15698,17523,18610,19842,
                20310,22500,23080,21916)
  
  r_less <- vonNeumannTest(tourists, alternative = "less")
  r_two  <- vonNeumannTest(tourists, alternative = "two.sided")
  
  expect_gte(r_two$p.value, r_less$p.value)
})

test_that("less and greater p-values sum to 1", {
  
  set.seed(1)
  x <- rnorm(20)
  
  r_less    <- vonNeumannTest(x, alternative = "less")
  r_greater <- vonNeumannTest(x, alternative = "greater")
  
  expect_equal(r_less$p.value + r_greater$p.value, 1, tolerance = 1e-10)
})

test_that("two-sided p = 2 * min(less, greater)", {
  
  set.seed(1)
  x <- rnorm(20)
  
  r_less    <- vonNeumannTest(x, alternative = "less")
  r_greater <- vonNeumannTest(x, alternative = "greater")
  r_two     <- vonNeumannTest(x, alternative = "two.sided")
  
  expect_equal(
    r_two$p.value,
    2 * min(r_less$p.value, r_greater$p.value),
    tolerance = 1e-10
  )
})

test_that("invalid alternative throws error", {
  expect_error(vonNeumannTest(rnorm(20), alternative = "invalid"))
})
# -------------------------------------------------------------------------
# NA handling
# -------------------------------------------------------------------------
test_that("NA values are silently removed", {
  
  set.seed(1)
  x    <- rnorm(25)
  x_na <- x
  x_na[c(3, 15)] <- NA
  
  res_clean <- vonNeumannTest(x[!is.na(x_na)])
  res_na    <- vonNeumannTest(x_na)
  
  expect_equal(unname(res_clean$statistic), unname(res_na$statistic))
  expect_equal(res_clean$p.value, res_na$p.value)
})
# -------------------------------------------------------------------------
# Input validation
# -------------------------------------------------------------------------
test_that("n < 4 throws error", {
  expect_error(vonNeumannTest(1:3), "4")
})

test_that("constant x throws error", {
  expect_error(vonNeumannTest(rep(1, 10)), "distinct")
})

test_that("non-numeric throws error", {
  expect_error(vonNeumannTest(letters[1:10]))
})

test_that("parameter n is correct", {
  
  set.seed(1)
  x <- rnorm(20)
  
  res <- vonNeumannTest(x)
  
  expect_equal(unname(res$parameter["n"]), 20L)
})

test_that("parameter n reflects NA removal", {
  
  set.seed(1)
  x    <- rnorm(22)
  x[c(1, 5)] <- NA
  
  res <- vonNeumannTest(x)
  
  expect_equal(unname(res$parameter["n"]), 20L)
})
# -------------------------------------------------------------------------
# Print compatibility
# -------------------------------------------------------------------------
test_that("print.htest works", {
  
  set.seed(2)
  x <- runif(20)
  
  res <- vonNeumannTest(x)
  
  expect_output(print(res), "Von Neumann")
})