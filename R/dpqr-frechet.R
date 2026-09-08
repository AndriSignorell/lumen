

#' Fréchet Distribution
#' 
#' The Fréchet distribution, also known as the Type II extreme value 
#' distribution, is a continuous probability distribution for the maximum 
#' of a sequence of independent random variables. It has a lower bound and 
#' a heavy right tail, and is parameterized by location, scale, and shape.
#' 
#' Density function, distribution function, quantile function and random
#' generation for the Frechet distribution with location, scale and shape
#' parameters.
#' 
#' The Frechet distribution function with parameters \eqn{`loc` = a},
#' \eqn{`scale` = b} and \eqn{`shape` = s} is 
#' \deqn{G(z) = \exp\left\{-\left(\frac{z-a}{b}\right)^{-s}\right\}}
#' for \eqn{z > a} and zero otherwise, where \eqn{b > 0} and \eqn{s > 0}.
#' 
#' @name dpqr-frechet
#' @aliases dfrechet pfrechet qfrechet rfrechet
#' 
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations.
#' @param loc,scale,shape location, scale and shape parameters (can be given as
#' vectors).
#' @param log,log.p logical; if `TRUE`, probabilities `p` are given as
#' `log(p)` and the density is returned on the log scale.
#' @param lower.tail logical; if `TRUE` (default), probabilities are 
#' \verb{P[X <= x]}, otherwise, P\verb{[X > x]}.
#' 
#' @return `dfrechet()` gives the density function, `pfrechet()`
#' gives the distribution function, `qfrechet()` gives the quantile
#' function, and `rfrechet()` generates random deviates.
#' 
#' @note
#' Based on code by Alec Stephenson previously published in
#' the \pkg{evd} package, adapted to conform to package standards.
#' 
#' @seealso [distributions-overview]
#' @concept distribution-function
#' @concept extreme-value
#' 
#' @examples
#' 
#' dfrechet(2:4, 1, 0.5, 0.8)
#' pfrechet(2:4, 1, 0.5, 0.8)
#' qfrechet(seq(0.9, 0.6, -0.1), 2, 0.5, 0.8)
#' rfrechet(6, 1, 0.5, 0.8)
#' p <- (1:9)/10
#' pfrechet(qfrechet(p, 1, 2, 0.8), 1, 2, 0.8)
#' ## [1] 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9
#' 


#' @rdname dpqr-frechet
#' @export
dfrechet <- function(x, loc = 0, scale = 1, shape = 1, log = FALSE)
  {
    .assertPositive(scale)
    .assertPositive(shape)
    x <- (x - loc)/scale
    xpos <- x[x>0 | is.na(x)]
    nn <- length(x)
    scale <- rep(scale, length.out = nn)[x>0 | is.na(x)]
    shape <- rep(shape, length.out = nn)[x>0 | is.na(x)]
    d <- numeric(nn)
    d[x>0 | is.na(x)] <- log(shape/scale) - (1+shape) * log(xpos) -
      xpos^(-shape)
    d[x<=0 & !is.na(x)] <- -Inf
    if(!log) d <- exp(d)
    d
  }

#' @rdname dpqr-frechet
#' @export
pfrechet <- function(q, loc = 0, scale = 1, shape = 1, lower.tail = TRUE,
                     log.p = FALSE)
  {
    .assertPositive(scale)
    .assertPositive(shape)
    q <- pmax((q - loc)/scale,0)
    p <- exp(-q^(-shape))
    if(!lower.tail) p <- 1 - p
    if(log.p) log(p) else p
  }

#' @rdname dpqr-frechet
#' @export
qfrechet <- function(p, loc = 0, scale = 1, shape = 1, lower.tail = TRUE,
                     log.p = FALSE)
  {
    .assertPositive(scale)
    .assertPositive(shape)
    p <- .qProb(p, lower.tail = lower.tail, log.p = log.p)
    loc + scale * (-log(p))^(-1/shape)
  }

#' @rdname dpqr-frechet
#' @export
rfrechet <- function(n, loc = 0, scale = 1, shape = 1)
  {
    .assertPositive(scale)
    .assertPositive(shape)
    loc + scale * rexp(n)^(-1/shape)
  }


