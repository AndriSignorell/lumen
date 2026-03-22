
#' The Generalized Extreme Value Distribution
#' 
#' The Generalized Extreme Value (GEV) distribution unifies the three 
#' extreme value distributions — Gumbel (Type I), Fréchet (Type II), 
#' and Reverse Weibull (Type III) — into a single family, parameterized 
#' by location, scale, and a shape parameter that determines which type applies.
#' 
#' Density function, distribution function, quantile function and random
#' generation for the generalized extreme value (GEV) distribution with
#' location, scale and shape parameters.
#' 
#' The GEV distribution function with parameters \eqn{\code{loc} = a},
#' \eqn{\code{scale} = b} and \eqn{\code{shape} = s} is \deqn{G(z) =
#' \exp\left[-\{1+s(z-a)/b\}^{-1/s}\right]}{ G(x) = exp[-{1+s(z-a)/b}^(-1/s)]}
#' for \eqn{1+s(z-a)/b > 0}, where \eqn{b > 0}.  If \eqn{s = 0} the
#' distribution is defined by continuity.  If \eqn{1+s(z-a)/b \leq
#' 0}{1+s(z-a)/b <= 0}, the value \eqn{z} is either greater than the upper end
#' point (if \eqn{s < 0}), or less than the lower end point (if \eqn{s > 0}).
#' 
#' The parametric form of the GEV encompasses that of the Gumbel, Frechet and
#' reverse Weibull distributions, which are obtained for \eqn{s = 0}, \eqn{s >
#' 0} and \eqn{s < 0} respectively.  It was first introduced by Jenkinson
#' (1955).
#' 
#' @name dgev
#' @aliases dgev pgev qgev rgev
#' @param x,q Vector of quantiles.
#' @param p Vector of probabilities.
#' @param n Number of observations.
#' @param loc,scale,shape Location, scale and shape parameters; the
#' \code{shape} argument cannot be a vector (must have length one).
#' @param log Logical; if \code{TRUE}, the log density is returned.
#' @param lower.tail Logical; if \code{TRUE} (default), probabilities are 
#' \verb{P[X <= x]}, otherwise, P\verb{[X > x]}
#' @return \code{dgev} gives the density function, \code{pgev} gives the
#' distribution function, \code{qgev} gives the quantile function, and
#' \code{rgev} generates random deviates.
#' 
#' @seealso \code{\link[evd]{fgev}}
#' 
#' @references Jenkinson, A. F. (1955) The frequency distribution of the annual
#' maximum (or minimum) of meteorological elements.  \emph{Quart. J. R. Met.
#' Soc.}, \bold{81}, 158--171.
#' 
#' @note
#' Based on code by Alec Stephenson. 
#' 
#' @family topic.extremevalue.distributions
#' @concept continuous distribution
#' @concept extreme value theory
#' @concept GEV
#' @concept dpqr
#' 
#' @examples
#' 
#' dgev(2:4, 1, 0.5, 0.8)
#' pgev(2:4, 1, 0.5, 0.8)
#' qgev(seq(0.9, 0.6, -0.1), 2, 0.5, 0.8)
#' rgev(6, 1, 0.5, 0.8)
#' p <- (1:9)/10
#' pgev(qgev(p, 1, 2, 0.8), 1, 2, 0.8)
#' ## [1] 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9
#' 

#' @rdname dgev
#' @export
"dgev"<-
  function(x, loc = 0, scale = 1, shape = 0, log = FALSE)
  {
    if(min(scale) <= 0) stop("invalid scale")
    if(length(shape) != 1) stop("invalid shape")
    x <- (x - loc)/scale
    if(shape == 0)
      d <- log(1/scale) - x - exp(-x) 
    else {
      nn <- length(x)
      xx <- 1 + shape*x
      xxpos <- xx[xx>0 | is.na(xx)]
      scale <- rep(scale, length.out = nn)[xx>0 | is.na(xx)]
      d <- numeric(nn)
      d[xx>0 | is.na(xx)] <- log(1/scale) - xxpos^(-1/shape) -
        (1/shape + 1)*log(xxpos)
      d[xx<=0 & !is.na(xx)] <- -Inf
    }  
    if(!log) d <- exp(d)
    d
  }


#' @rdname dgev
#' @export
"pgev"<-
  function(q, loc = 0, scale = 1, shape = 0, lower.tail = TRUE)
  {
    if(min(scale) <= 0) stop("invalid scale")
    if(length(shape) != 1) stop("invalid shape")
    q <- (q - loc)/scale
    if(shape == 0) p <- exp(-exp(-q))
    else p <- exp( - pmax(1 + shape * q, 0)^(-1/shape))
    if(!lower.tail) p <- 1 - p
    p
  }


#' @rdname dgev
#' @export
"qgev"<-
  function(p, loc = 0, scale = 1, shape = 0, lower.tail = TRUE)
  {
    if(min(p, na.rm = TRUE) <= 0 || max(p, na.rm = TRUE) >=1)
      stop("`p' must contain probabilities in (0,1)")
    if(min(scale) < 0) stop("invalid scale")
    if(length(shape) != 1) stop("invalid shape")
    if(!lower.tail) p <- 1 - p
    if(shape == 0) return(loc - scale * log(-log(p)))
    else return(loc + scale * ((-log(p))^(-shape) - 1)/shape)
  }


#' @rdname dgev
#' @export
"rgev"<-
  function(n, loc = 0, scale = 1, shape = 0)
  {
    if(min(scale) < 0) stop("invalid scale")
    if(length(shape) != 1) stop("invalid shape")
    if(shape == 0) return(loc - scale * log(rexp(n)))
    else return(loc + scale * (rexp(n)^(-shape) - 1)/shape)
  }



