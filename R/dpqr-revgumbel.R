
#' Reverse Gumbel Distribution
#' 
#' Density, distribution function, quantile function and random generation for
#' the \dQuote{Reverse} Gumbel distribution with parameters `loc` and
#' `scale`.
#' 
#' The reverse Gumbel distribution is the distribution of \eqn{a - bY} for a
#' standard Gumbel \eqn{Y}, i.e. the Type I extreme value distribution for
#' minima. With \eqn{`loc` = a} and \eqn{`scale` = b} its distribution
#' function is
#' \deqn{F(x) = 1 - \exp\left\{-\exp\left[\left(\frac{x-a}{b}\right)\right]\right\}}{F(x) = 1 - exp(-exp((x-a)/b))}
#' for all real \eqn{x}, where \eqn{b > 0}.
#' 
#' @name dpqr-revgumbel
#' @aliases drevgumbel prevgumbel qrevgumbel qrevgumbelExp rrevgumbel
#' @param x,q numeric vector of abscissa (or quantile) values at which to
#' evaluate the density or distribution function.
#' @param p numeric vector of probabilities at which to evaluate the quantile
#' function.
#' @param loc location of the distribution.
#' @param scale scale (\eqn{> 0}) of the distribution.
#' @param n number of random variates, i.e., [length()] of resulting
#' vector of `rrevgumbel()`.
#' @param log,log.p logical; if `TRUE`, probabilities `p` are given as
#' `log(p)` and the density is returned on the log scale.
#' @param lower.tail logical; if `TRUE` (default), probabilities are 
#' \verb{P[X <= x]}, otherwise, P\verb{[X > x]}.
#' @return A numeric vector, of the same length as `x`, `q`, or
#' `p` for the first three functions, and of length `n` for
#' `rrevgumbel()`. `qrevgumbelExp()` gives the quantiles of
#' \eqn{\exp(X)}, the exponential parametrization used in some applications.
#' @seealso [distributions-overview]; [dpqr-gumbel] for the Gumbel
#' distribution this one reverses.
#' @note
#' Based on code by Werner Stahel, partly inspired by the \pkg{VGAM} package
#' (numeric refinements by Martin Maechler), adapted to conform to package
#' standards.
#' 
#' @examples
#' 
#' curve(prevgumbel(x, scale= 1/2), -3,2, n=1001, col=1, lwd=2,
#'       main = "revgumbel(x, scale = 1/2)")
#' abline(h=0:1, v = 0, lty=3, col = "gray30")
#' curve(drevgumbel(x, scale= 1/2),       n=1001, add=TRUE,
#'       col = (col.d <- adjustcolor(2, 0.5)), lwd=3)
#' legend("left", c("cdf","pdf"), col=c("black", col.d), lwd=2:3, bty="n")
#' 
#' med <- qrevgumbel(0.5, scale=1/2)
#' cat("The median is:",  format(med),"\n")
#' 

#' @rdname dpqr-revgumbel
#' @concept distribution-function
#' @concept extreme-value
#' @export
drevgumbel <- function (x, loc = 0, scale = 1, log = FALSE) {
  .assertPositive(scale)
  t <- (x - loc)/scale
  d <- t - exp(t) - log(scale)
  if (log) d else exp(d)
}

#' @rdname dpqr-revgumbel
#' @export
prevgumbel <- function (q, loc = 0, scale = 1, lower.tail = TRUE,
                        log.p = FALSE) {
  .assertPositive(scale)
  # the upper tail is exp(-exp(t)) and thus exact on the log scale
  lupper <- -exp((q - loc)/scale)
  if (lower.tail) {
    if (log.p) log(-expm1(lupper)) else -expm1(lupper)
  } else {
    if (log.p) lupper else exp(lupper)
  }
}

#' @rdname dpqr-revgumbel
#' @export
qrevgumbel <- function (p, loc = 0, scale = 1, lower.tail = TRUE,
                        log.p = FALSE) {
  .assertPositive(scale)
  p <- .qProb(p, lower.tail = lower.tail, log.p = log.p)
  loc + scale * log(-log1p(-p))
}

#' @rdname dpqr-revgumbel
#' @export
rrevgumbel <- function (n, loc = 0, scale = 1) {
  .assertPositive(scale)
  loc + scale * log(-log(runif(n)))
}

#' @rdname dpqr-revgumbel
#' @export
qrevgumbelExp <- function (p, loc = 0, scale = 1, lower.tail = TRUE,
                           log.p = FALSE)
  exp(qrevgumbel(p, loc = loc, scale = scale, lower.tail = lower.tail,
                 log.p = log.p))
