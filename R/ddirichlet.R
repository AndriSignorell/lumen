#' Dirichlet distribution functions
#'
#' Density, random generation, and basic utilities for the Dirichlet distribution.
#'
#' The Dirichlet distribution is a multivariate generalization of the Beta
#' distribution defined on the simplex:
#' \deqn{\sum_{i=1}^k x_i = 1, \quad x_i \ge 0}
#'
#' @name dirichlet
#' 
#' @family multivariate-distributions
#' @family dirichlet
#' @concept multivariate-distribution
#' @concept simplex
#' @concept dirichlet-distribution
#' @concept bayesian
#' @concept dpqr
NULL


#' Dirichlet density
#'
#' Computes the density of the Dirichlet distribution.
#'
#' @param x Numeric vector or matrix (rows sum to 1)
#' @param alpha Numeric vector of concentration parameters (> 0)
#' @param log Logical; return log-density if TRUE
#'
#' @return Numeric vector of densities
#'
#' @examples
#' ddirichlet(c(0.2, 0.3, 0.5), c(1,1,1))
#'
#' @family dist.other
#' @concept distributions
#' @concept multivariate
#'
#'
#' @export
ddirichlet <- function(x, alpha, log = FALSE) {
  if (is.vector(x)) x <- matrix(x, nrow = 1)
  
  if (any(alpha <= 0)) {
    stop("alpha must be > 0")
  }
  
  if (ncol(x) != length(alpha)) {
    stop("x and alpha must have same length")
  }
  
  # check simplex (soft check)
  if (any(x < 0) || any(abs(rowSums(x) - 1) > 1e-8)) {
    return(rep(if (log) -Inf else 0, nrow(x)))
  }
  
  logdens <- lgamma(sum(alpha)) - sum(lgamma(alpha)) +
    rowSums((alpha - 1) * log(x))
  
  if (log) logdens else exp(logdens)
}



#' Dirichlet CDF (parallel Monte Carlo)
#'
#' Fast parallel approximation of the Dirichlet CDF using RcppParallel.
#'
#' @param q Numeric vector
#' @param alpha Numeric vector (> 0)
#' @param n_sim Number of simulations
#'
#' @return Approximate probability
#'
#' @examples
#' pdirichlet(c(0.2,0.3,0.5), c(1,1,1))
#'

#' @export
pdirichlet <- function(q, alpha, n_sim = 1e5) {
  pdirichlet_cpp(q, alpha, as.integer(n_sim))
}


#' Random Dirichlet
#'
#' Generates random draws from a Dirichlet distribution.
#'
#' @param n Number of samples
#' @param alpha Numeric vector of concentration parameters (> 0)
#'
#' @return Matrix with n rows
#'
#' @examples
#' rdirichlet(5, c(1,1,1))
#'
#' @export
rdirichlet <- function(n, alpha) {
  if (any(alpha <= 0)) {
    stop("alpha must be > 0")
  }
  
  k <- length(alpha)
  
  x <- matrix(rgamma(n * k, shape = alpha, rate = 1), nrow = n)
  x / rowSums(x)
}



#' Dirichlet quantile (not defined)
#'
#' The Dirichlet distribution has no unique multivariate quantile function.
#'
#' @export
qdirichlet <- function() {
  stop("Quantile function is not defined for the Dirichlet distribution.")
}


