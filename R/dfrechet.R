

#' The Frechet Distribution
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
#' The Frechet distribution function with parameters \eqn{\code{loc} = a},
#' \eqn{\code{scale} = b} and \eqn{\code{shape} = s} is 
#' \deqn{G(z) = \exp\left\{-\left(\frac{z-a}{b}\right)^{-s}\right\}}
#' for \eqn{z > a} and zero otherwise, where \eqn{b > 0} and \eqn{s > 0}.
#' 
#' @name dfrechet
#' @aliases dfrechet pfrechet qfrechet rfrechet
#' @param x,q Vector of quantiles.
#' @param p Vector of probabilities.
#' @param n Number of observations.
#' @param loc,scale,shape Location, scale and shape parameters (can be given as
#' vectors).
#' @param log Logical; if \code{TRUE}, the log density is returned.
#' @param lower.tail Logical; if \code{TRUE} (default), probabilities are 
#' \verb{P[X <= x]}, otherwise, P\verb{[X > x]}
#' @return \code{dfrechet} gives the density function, \code{pfrechet} gives
#' the distribution function, \code{qfrechet} gives the quantile function, and
#' \code{rfrechet} generates random deviates.
#' 
#' @note
#' Based on code by Alec Stephenson. 
#' 
#' @family topic.extremevalue.distributions
#' @concept continuous distribution
#' @concept extreme value theory
#' @concept Frechet
#' @concept heavy tail
#' @concept dpqr
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


#' @rdname dfrechet
#' @export
"dfrechet"<-
  function(x, loc = 0, scale = 1, shape = 1, log = FALSE)
  {
    if(min(scale) <= 0 || min(shape) <= 0) stop("invalid arguments")
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

#' @rdname dfrechet
#' @export
"pfrechet"<-
  function(q, loc = 0, scale = 1, shape = 1, lower.tail = TRUE)
  {
    if(min(scale) <= 0 || min(shape) <= 0) stop("invalid arguments")
    q <- pmax((q - loc)/scale,0)
    p <- exp(-q^(-shape))
    if(!lower.tail) p <- 1 - p
    p
  }

#' @rdname dfrechet
#' @export
"qfrechet"<-
  function(p, loc = 0, scale = 1, shape = 1, lower.tail = TRUE)
  {
    if(min(p, na.rm = TRUE) <= 0 || max(p, na.rm = TRUE) >=1)
      stop("`p' must contain probabilities in (0,1)")
    if(min(scale) < 0 || min(shape) <= 0) stop("invalid arguments")
    if(!lower.tail) p <- 1 - p
    loc + scale * (-log(p))^(-1/shape)
  }

#' @rdname dfrechet
#' @export
"rfrechet"<-
  function(n, loc = 0, scale = 1, shape = 1)
  {
    if(min(scale) < 0 || min(shape) <= 0) stop("invalid arguments")
    loc + scale * rexp(n)^(-1/shape)
  }


