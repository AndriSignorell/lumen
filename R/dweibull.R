
#' The Reverse Weibull Distribution
#' 
#' Density function, distribution function, quantile function and random
#' generation for the reverse Weibull and the negative Weibull distribution with 
#' location, scale and shape parameters.
#' 
#' The reverse (or negative) Weibull distribution function with parameters
#' \eqn{\code{loc} = a}, \eqn{\code{scale} = b} and \eqn{\code{shape} = s} is
#' \deqn{G(z) = \exp\left\{-\left[-\left(\frac{z-a}{b}\right)\right]^s\right\}}{G(z) = exp(-(-(z-a)/b)^s)}
#' for \eqn{z < a} and one otherwise, where \eqn{b > 0} and \eqn{s > 0}.
#'  
#' @name drweibull
#' @aliases drweibull prweibull qrweibull rrweibull dnweibull pnweibull
#' qnweibull rnweibull
#' @param x,q Vector of quantiles.
#' @param p Vector of probabilities.
#' @param n Number of observations.
#' @param loc,scale,shape Location, scale and shape parameters (can be given as
#' vectors).
#' @param log Logical; if \code{TRUE}, the log density is returned.
#' @param lower.tail Logical; if \code{TRUE} (default), probabilities are 
#' \verb{P[X <= x]}, otherwise, P\verb{[X > x]}
#' @return \code{drweibull} and \code{dnweibull} give the density function,
#' \code{prweibull} and \code{pnweibull} give the distribution function,
#' \code{qrweibull} and \code{qnweibull} give the quantile function,
#' \code{rrweibull} and \code{rnweibull} generate random deviates.
#' @note Within extreme value theory the reverse Weibull distibution (also
#' known as the negative Weibull distribution) is often referred to as the
#' Weibull distribution.  We make a distinction to avoid confusion with the
#' three-parameter distribution used in survival analysis, which is related by
#' a change of sign to the distribution given above.
#' @seealso \code{\link{rfrechet}}, \code{\link{rgev}}, \code{\link{rgumbel}}
#' 
#' @family topic.distributions
#' @concept continuous distribution
#' @concept extreme value theory
#' @concept reverse Weibull
#' @concept bounded tail
#' @concept dpqr
#' 
#' @examples
#' 
#' drweibull(-5:-3, -1, 0.5, 0.8)
#' prweibull(-5:-3, -1, 0.5, 0.8)
#' qrweibull(seq(0.9, 0.6, -0.1), 2, 0.5, 0.8)
#' rrweibull(6, -1, 0.5, 0.8)
#' p <- (1:9)/10
#' prweibull(qrweibull(p, -1, 2, 0.8), -1, 2, 0.8)
#' ## [1] 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9
#' 


#' @rdname drweibull
#' @export
"drweibull"<-
  function(x, loc = 0, scale = 1, shape = 1, log = FALSE)
  {
    if(min(scale) <= 0 || min(shape) <= 0) stop("invalid arguments")
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


#' @rdname drweibull
#' @export
"prweibull"<-
  function(q, loc = 0, scale = 1, shape = 1, lower.tail = TRUE)
  {
    if(min(scale) <= 0 || min(shape) <= 0) stop("invalid arguments")
    q <- pmin((q - loc)/scale,0)
    p <- exp(-(-q)^shape)
    if(!lower.tail) p <- 1 - p
    p
  }

#' @rdname drweibull
#' @export
"qrweibull"<-
  function(p, loc = 0, scale = 1, shape = 1, lower.tail = TRUE)
  {
    if(min(p, na.rm = TRUE) <= 0 || max(p, na.rm = TRUE) >=1)
      stop("`p' must contain probabilities in (0,1)")
    if(min(scale) < 0 || min(shape) <= 0) stop("invalid arguments")
    if(!lower.tail) p <- 1 - p
    loc - scale * (-log(p))^(1/shape)
  }

#' @rdname drweibull
#' @export
"rrweibull"<-
  function(n, loc = 0, scale = 1, shape = 1)
  {
    if(min(scale) < 0 || min(shape) <= 0) stop("invalid arguments")
    loc - scale * rexp(n)^(1/shape)
  }





"dnweibull"<-
  function(x, loc = 0, scale = 1, shape = 1, log = FALSE)
  {
    if(min(scale) <= 0 || min(shape) <= 0) stop("invalid arguments")
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

"pnweibull"<-
  function(q, loc = 0, scale = 1, shape = 1, lower.tail = TRUE)
  {
    if(min(scale) <= 0 || min(shape) <= 0) stop("invalid arguments")
    q <- pmin((q - loc)/scale,0)
    p <- exp(-(-q)^shape)
    if(!lower.tail) p <- 1 - p
    p
  }


"qnweibull"<-
  function(p, loc = 0, scale = 1, shape = 1, lower.tail = TRUE)
  {
    if(min(p, na.rm = TRUE) <= 0 || max(p, na.rm = TRUE) >=1)
      stop("`p' must contain probabilities in (0,1)")
    if(min(scale) < 0 || min(shape) <= 0) stop("invalid arguments")
    if(!lower.tail) p <- 1 - p
    loc - scale * (-log(p))^(1/shape)
  }

"rnweibull"<-
  function(n, loc = 0, scale = 1, shape = 1)
  {
    if(min(scale) < 0 || min(shape) <= 0) stop("invalid arguments")
    loc - scale * rexp(n)^(1/shape)
  }

