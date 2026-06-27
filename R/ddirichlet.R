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
#' @param concentration Numeric vector of concentration parameters (> 0)
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
ddirichlet <- function(x, concentration, log = FALSE) {
  if (is.vector(x)) x <- matrix(x, nrow = 1)
  
  if (any(concentration <= 0)) {
    stop("concentration must be > 0")
  }
  
  if (ncol(x) != length(concentration)) {
    stop("x and concentration must have same length")
  }
  
  # check simplex (soft check)
  if (any(x < 0) || any(abs(rowSums(x) - 1) > 1e-8)) {
    return(rep(if (log) -Inf else 0, nrow(x)))
  }
  
  logdens <- lgamma(sum(concentration)) - sum(lgamma(concentration)) +
    rowSums((concentration - 1) * log(x))
  
  if (log) logdens else exp(logdens)
}



#' Dirichlet CDF (parallel Monte Carlo)
#'
#' Fast parallel approximation of the Dirichlet CDF using RcppParallel.
#'
#' @param q Numeric vector
#' @param concentration Numeric vector (> 0)
#' @param nSim Number of simulations
#'
#' @return Approximate probability
#'
#' @examples
#' pdirichlet(c(0.2,0.3,0.5), c(1,1,1))
#'

#' @export
pdirichlet <- function(q, concentration, nSim = 1e5) {
  pdirichlet_cpp(q, concentration, as.integer(nSim))
}


#' Random Dirichlet
#'
#' Generates random draws from a Dirichlet distribution.
#'
#' @param n Number of samples
#' @param concentration Numeric vector of concentration parameters (> 0)
#'
#' @return Matrix with n rows
#'
#' @examples
#' rdirichlet(5, c(1,1,1))
#'
#' @export
rdirichlet <- function(n, concentration) {
  if (any(concentration <= 0)) {
    stop("concentration must be > 0")
  }
  
  k <- length(concentration)
  
  x <- matrix(rgamma(n * k, shape = concentration, rate = 1), nrow = n)
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


