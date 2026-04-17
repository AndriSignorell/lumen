
library(testthat)

# ------------------------------------------------------------------------------
# basic functionality
# ------------------------------------------------------------------------------

test_that("pdirichlet returns scalar in [0,1]", {
  q <- c(0.2, 0.3, 0.5)
  alpha <- c(1, 1, 1)
  
  p <- pdirichlet(q, alpha, n_sim = 1e4)
  
  expect_length(p, 1)
  expect_true(is.numeric(p))
  expect_true(p >= 0 && p <= 1)
})

# ------------------------------------------------------------------------------
# symmetry (Dirichlet(1,1,...,1))
# ------------------------------------------------------------------------------

test_that("symmetry for uniform Dirichlet", {
  q1 <- c(0.2, 0.3, 0.5)
  q2 <- c(0.5, 0.3, 0.2)
  
  alpha <- c(1, 1, 1)
  
  p1 <- pdirichlet(q1, alpha, n_sim = 2e4)
  p2 <- pdirichlet(q2, alpha, n_sim = 2e4)
  
  expect_equal(p1, p2, tolerance = 0.02)
})

# ------------------------------------------------------------------------------
# monotonicity
# ------------------------------------------------------------------------------

test_that("CDF is monotone in q", {
  alpha <- c(1,1,1)
  
  q_small <- c(0.1, 0.2, 0.7)
  q_large <- c(0.2, 0.3, 0.5)
  
  p_small <- pdirichlet(q_small, alpha, n_sim = 2e4)
  p_large <- pdirichlet(q_large, alpha, n_sim = 2e4)
  
  expect_true(p_small <= p_large)
})

# ------------------------------------------------------------------------------
# extreme cases
# ------------------------------------------------------------------------------

test_that("extreme q values", {
  alpha <- c(1,1,1)
  
  # very small region
  q_small <- c(0.01, 0.01, 0.98)
  p_small <- pdirichlet(q_small, alpha, n_sim = 2e4)
  
  expect_true(p_small < 0.2)
  
  # large region (almost everything)
  q_large <- c(0.99, 0.99, 0.99)
  p_large <- pdirichlet(q_large, alpha, n_sim = 2e4)
  
  expect_true(p_large > 0.95)
})

# ------------------------------------------------------------------------------
# dimension handling
# ------------------------------------------------------------------------------

test_that("works for higher dimensions", {
  alpha <- rep(1, 5)
  q <- rep(1/5, 5)
  
  p <- pdirichlet(q, alpha, n_sim = 2e4)
  
  expect_true(p >= 0 && p <= 1)
})

# ------------------------------------------------------------------------------
# invalid input
# ------------------------------------------------------------------------------

test_that("invalid input throws error", {
  expect_error(
    pdirichlet(c(0.2,0.8), c(1,1,1), n_sim = 1000)
  )
})

# ------------------------------------------------------------------------------
# consistency with slow R version
# ------------------------------------------------------------------------------

test_that("agrees with R Monte Carlo implementation", {
  
  pdirichlet_R <- function(q, alpha, n_sim = 1e4) {
    sims <- rdirichlet(n_sim, alpha)
    mean(apply(sims <= matrix(q, nrow = n_sim, ncol = length(q), byrow = TRUE), 1, all))
  }
  
  q <- c(0.2, 0.3, 0.5)
  alpha <- c(1,1,1)
  
  p_cpp <- pdirichlet(q, alpha, n_sim = 2e4)
  p_r   <- pdirichlet_R(q, alpha, n_sim = 2e4)
  
  expect_equal(p_cpp, p_r, tolerance = 0.02)
})