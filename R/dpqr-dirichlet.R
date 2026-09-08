#' Dirichlet Distribution
#'
#' Density, random generation, and basic utilities for the Dirichlet distribution.
#'
#' The Dirichlet distribution is a multivariate generalization of the Beta
#' distribution defined on the simplex:
#' \deqn{\sum_{i=1}^k x_i = 1, \quad x_i \ge 0}
#'
#' @name dpqr-dirichlet
#' @aliases ddirichlet pdirichlet qdirichlet rdirichlet
#' 
#' @param x numeric vector or matrix (rows sum to 1).
#' @param q numeric vector of quantiles.
#' @param p numeric vector of probabilities; accepted by `qdirichlet()`
#' only to give a meaningful error, see below.
#' @param n number of samples.
#' @param concentration numeric vector of concentration parameters (> 0).
#' @param R number of Monte Carlo simulations used to approximate the CDF.
#'   Must be at most `.Machine$integer.max`.
#' @param log logical; return log-density if TRUE.
#' @param \dots further arguments, accepted by `qdirichlet()` and ignored.
#'
#' @section Random number generation:
#' `pdirichlet()` evaluates the CDF by simulation and is parallelised. Each
#' range of draws runs its own generator, seeded from R's stream, so that
#' [set.seed()] governs the result; the seeding does depend on how the work
#' is split, so reproducing a value also requires the same
#' `RcppParallel::setThreadOptions()`.
#'
#' @return `ddirichlet()` gives a numeric vector of densities (one per
#' row of `x`), `pdirichlet()` gives an approximate probability,
#' and `rdirichlet()` generates a matrix with `n` rows of random
#' deviates. `qdirichlet()` only signals an error, as no unique
#' multivariate quantile function exists; it takes the same arguments as
#' the others so that the message is reached rather than an argument
#' mismatch.
#'
#' @examples
#' ddirichlet(c(0.2, 0.3, 0.5), c(1,1,1))
#' pdirichlet(c(0.2, 0.3, 0.5), c(1,1,1))
#' rdirichlet(5, c(1,1,1))
#'
#' @seealso [distributions-overview]
#' @concept distribution-function
#' @concept multivariate
#' @concept proportion
#'
#' @rdname dpqr-dirichlet
#' @export
ddirichlet <- function(x, concentration, log = FALSE) {
  if (is.vector(x)) x <- matrix(x, nrow = 1)

  .assertPositive(concentration)

  if (ncol(x) != length(concentration))
    stop("'x' and 'concentration' must have the same length", call. = FALSE)
  
  w <- concentration - 1
  lx <- log(pmax(x, 0))
  lx[, w == 0] <- 0   # avoid 0 * log(0) = NaN on the simplex boundary
  
  logdens <- as.vector(lx %*% w) +
    lgamma(sum(concentration)) - sum(lgamma(concentration))
  
  # rows off the simplex have density 0 (soft check, rowwise)
  rs <- rowSums(x)
  offSimplex <- apply(x, 1, function(r) isTRUE(any(r < 0))) |
    (!is.na(rs) & abs(rs - 1) > 1e-8)
  logdens[offSimplex] <- -Inf
  
  if (log) logdens else exp(logdens)
}



#' @rdname dpqr-dirichlet
#' @export
pdirichlet <- function(q, concentration, R = 1e5) {
  .assertPositive(concentration)
  .assertScalar(R, lower = 1, upper = .Machine$integer.max,
                integerValued = TRUE)

  if (length(q) != length(concentration))
    stop("'q' and 'concentration' must have the same length", call. = FALSE)

  pdirichlet_cpp(q, concentration, as.integer(R))
}


#' @rdname dpqr-dirichlet
#' @export
rdirichlet <- function(n, concentration) {
  .assertPositive(concentration)

  k <- length(concentration)
  
  # byrow = TRUE so that each row receives shapes concentration[1..k]
  # (rgamma recycles 'shape' sequentially over the n*k draws)
  x <- matrix(rgamma(n * k, shape = concentration, rate = 1),
              ncol = k, byrow = TRUE)
  x / rowSums(x)
}



#' @rdname dpqr-dirichlet
#' @export
qdirichlet <- function(p, concentration, ...) {
  stop("no quantile function is defined for the Dirichlet distribution",
       call. = FALSE)
}


