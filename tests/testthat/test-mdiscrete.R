library(testthat)
library(lumen)

tol <- 1e-10

# --- mbinom ---
test_that("mbinom: mean=n*p, variance=n*p*(1-p)", {
  r <- mbinom(10, 0.3)
  expect_equal(unname(r["mean"]),     3,   tolerance = tol)
  expect_equal(unname(r["variance"]), 2.1, tolerance = tol)
})

# --- mpois ---
test_that("mpois: mean=variance=lambda", {
  r <- mpois(4)
  expect_equal(unname(r["mean"]),     4, tolerance = tol)
  expect_equal(unname(r["variance"]), 4, tolerance = tol)
})

# --- mgeom ---
test_that("mgeom: mean=(1-p)/p, variance=(1-p)/p^2", {
  r <- mgeom(0.5)
  expect_equal(unname(r["mean"]),     1,   tolerance = tol)
  expect_equal(unname(r["variance"]), 2,   tolerance = tol)
})

# --- mnbinom ---
test_that("mnbinom: mean=size*(1-p)/p, variance=mean/p", {
  r <- mnbinom(5, 0.5)
  expect_equal(unname(r["mean"]),     5,  tolerance = tol)
  expect_equal(unname(r["variance"]), 10, tolerance = tol)
})

# --- mhyper ---
test_that("mhyper: mean=k*m/(m+n)", {
  # mhyper(m, n, k): m black, n white, k drawn
  r <- mhyper(10, 10, 5)
  expect_equal(unname(r["mean"]), 2.5, tolerance = tol)
})

# --- mbenford ---
test_that("mbenford: mean and variance are numeric and positive", {
  r <- mbenford(1)
  expect_true(is.numeric(r["mean"]) && r["mean"] > 0)
  expect_true(is.numeric(r["variance"]) && r["variance"] > 0)
})

test_that("mbenford: ndigits=2 gives mean in [10,99]", {
  r <- mbenford(2)
  expect_gt(unname(r["mean"]), 10)
  expect_lt(unname(r["mean"]), 99)
})

# --- return structure ---
test_that("all discrete moment functions return mean and variance", {
  fns <- list(mbinom(10,0.5), mpois(3), mgeom(0.4), mnbinom(3, 0.5))
  for (f in fns) {
    expect_true(all(c("mean","variance") %in% names(f)))
  }
})
