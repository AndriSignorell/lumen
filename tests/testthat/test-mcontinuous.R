library(testthat)
library(lumen)

tol <- 1e-10

# --- mnorm ---
test_that("mnorm: mean and variance correct", {
  r <- mnorm(3, 2)
  expect_equal(unname(r["mean"]),     3,   tolerance = tol)
  expect_equal(unname(r["variance"]), 4,   tolerance = tol)
})

# --- mexp ---
test_that("mexp: mean=1/rate, variance=1/rate^2", {
  r <- mexp(2)
  expect_equal(unname(r["mean"]),     0.5,  tolerance = tol)
  expect_equal(unname(r["variance"]), 0.25, tolerance = tol)
})

# --- mgamma ---
test_that("mgamma: mean=shape/rate, variance=shape/rate^2", {
  r <- mgamma(4, 2)
  expect_equal(unname(r["mean"]),     2,   tolerance = tol)
  expect_equal(unname(r["variance"]), 1,   tolerance = tol)
})

# --- mlnorm ---
test_that("mlnorm: mean=exp(mu+sig^2/2), variance formula", {
  r <- mlnorm(0, 1)
  expect_equal(unname(r["mean"]),     exp(0.5),           tolerance = tol)
  expect_equal(unname(r["variance"]), (exp(1)-1)*exp(1),  tolerance = tol)
})

# --- mbeta ---
test_that("mbeta: mean=a/(a+b), variance formula", {
  r <- mbeta(2, 3)
  expect_equal(unname(r["mean"]),     2/5,                    tolerance = tol)
  expect_equal(unname(r["variance"]), 2*3/(25*6),             tolerance = tol)
})

# --- mchisq ---
test_that("mchisq: mean=df, variance=2*df", {
  r <- mchisq(5)
  expect_equal(unname(r["mean"]),     5,  tolerance = tol)
  expect_equal(unname(r["variance"]), 10, tolerance = tol)
})

# --- mt ---
test_that("mt: mean=0, variance=df/(df-2) for df>2", {
  r <- mt(10)
  expect_equal(unname(r["mean"]),     0,         tolerance = tol)
  expect_equal(unname(r["variance"]), 10/8,      tolerance = tol)
})

# --- mf ---
test_that("mf: mean=df2/(df2-2) for df2>2", {
  r <- mf(5, 10)
  expect_equal(unname(r["mean"]),     10/8, tolerance = tol)
})

# --- mtri ---
test_that("mtri: mean=(min+max+mode)/3", {
  r <- mtri(0, 6, 3)
  expect_equal(unname(r["mean"]),     3,   tolerance = tol)
  expect_equal(unname(r["variance"]), 3/2, tolerance = tol)
})

# --- return structure ---
test_that("all continuous moment functions return mean and variance", {
  fns <- list(mnorm(0,1), mexp(1), mgamma(2,1), mchisq(4), mt(5), mf(4,10))
  for (f in fns) {
    expect_true(all(c("mean","variance") %in% names(f)))
  }
})
