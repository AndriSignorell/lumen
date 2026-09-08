
#' Reverse Weibull Distribution
#' 
#' The Reverse Weibull distribution, also known as the Type III extreme 
#' value distribution, is the distribution of the negative of a 
#' Weibull-distributed random variable. It has an upper bound and a 
#' left-skewed density, and is parameterized by location, scale, and shape.
#' 
#' Density function, distribution function, quantile function and random
#' generation for the reverse (sometimes called negative) Weibull distribution with 
#' location, scale and shape parameters.
#' 
#' The reverse Weibull distribution function with parameters
#' \eqn{`loc` = a}, \eqn{`scale` = b} and \eqn{`shape` = s} is
#' \deqn{G(z) = \exp\left\{-\left[-\left(\frac{z-a}{b}\right)\right]^s\right\}}{G(z) = exp(-(-(z-a)/b)^s)}
#' for \eqn{z < a} and one otherwise, where \eqn{b > 0} and \eqn{s > 0}.
#'  
#' **Note:** Within extreme value theory the reverse Weibull distibution (also
#' known as the negative Weibull distribution) is often referred to as the
#' Weibull distribution.  We make a distinction to avoid confusion with the
#' three-parameter distribution used in survival analysis, which is related by
#' a change of sign to the distribution given above.
#' 
#' @name dpqr-revweibull
#' @aliases drevweibull prevweibull qrevweibull rrevweibull dnweibull pnweibull qnweibull rnweibull
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
#' @return `drevweibull()` and `dnweibull()` give the density
#' function, `prevweibull()` and `pnweibull()` give the distribution
#' function, `qrevweibull()` and `qnweibull()` give the quantile
#' function, `rrevweibull()` and `rnweibull()` generate random
#' deviates.
#' @seealso [distributions-overview]
#' 
#' @note
#' Based on code by Alec Stephenson previously published in
#' the \pkg{evd} package, adapted to conform to package standards.
#' 
#' @concept distribution-function
#' @concept extreme-value
#' 
#' @examples
#' 
#' drevweibull(-5:-3, -1, 0.5, 0.8)
#' prevweibull(-5:-3, -1, 0.5, 0.8)
#' qrevweibull(seq(0.9, 0.6, -0.1), 2, 0.5, 0.8)
#' rrevweibull(6, -1, 0.5, 0.8)
#' p <- (1:9)/10
#' prevweibull(qrevweibull(p, -1, 2, 0.8), -1, 2, 0.8)
#' ## [1] 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9
#' 


#' @rdname dpqr-revweibull
#' @export
drevweibull <- function(x, loc = 0, scale = 1, shape = 1, log = FALSE)
  {
    .assertPositive(scale)
    .assertPositive(shape)
    x <- (x - loc)/scale
    xneg <- x[x<0 | is.na(x)]
    nn <- length(x)
    scale <- rep(scale, length.out = nn)[x<0 | is.na(x)]
    shape <- rep(shape, length.out = nn)[x<0 | is.na(x)]
    d <- numeric(nn)
    d[x<0 | is.na(x)] <- log(shape/scale) + (shape-1) * log(-xneg) -
      (-xneg)^shape
    d[x>=0 & !is.na(x)] <- -Inf
    if(!log) d <- exp(d)
    d
  }



#' @rdname dpqr-revweibull
#' @export
prevweibull <- function(q, loc = 0, scale = 1, shape = 1, lower.tail = TRUE,
                        log.p = FALSE)
  {
    .assertPositive(scale)
    .assertPositive(shape)
    q <- pmin((q - loc)/scale,0)
    p <- exp(-(-q)^shape)
    if(!lower.tail) p <- 1 - p
    if(log.p) log(p) else p
  }

#' @rdname dpqr-revweibull
#' @export
qrevweibull <- function(p, loc = 0, scale = 1, shape = 1, lower.tail = TRUE,
                        log.p = FALSE)
  {
    .assertPositive(scale)
    .assertPositive(shape)
    p <- .qProb(p, lower.tail = lower.tail, log.p = log.p)
    loc - scale * (-log(p))^(1/shape)
  }

#' @rdname dpqr-revweibull
#' @export
rrevweibull <- function(n, loc = 0, scale = 1, shape = 1)
  {
    .assertPositive(scale)
    .assertPositive(shape)
    loc - scale * rexp(n)^(1/shape)
  }


# "negative Weibull" is an alternative name for the same reverse Weibull
# distribution (see Details) -- exported as plain synonyms, not separate
# implementations.

#' @rdname dpqr-revweibull
#' @export
dnweibull <- drevweibull

#' @rdname dpqr-revweibull
#' @export
pnweibull <- prevweibull

#' @rdname dpqr-revweibull
#' @export
qnweibull <- qrevweibull

#' @rdname dpqr-revweibull
#' @export
rnweibull <- rrevweibull



