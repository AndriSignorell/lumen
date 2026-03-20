
#' Maxima of Two Gumbel Distributions
#' 
#' Density function, distribution function, quantile function and random
#' generation for the maxima of two Gumbel distributions, each with different
#' location and scale parameters.
#' 
#' 
#' @name dgumbelx
#' @aliases dgumbelx pgumbelx qgumbelx rgumbelx
#' @param x,q Vector of quantiles.
#' @param p Vector of probabilities.
#' @param n Number of observations.
#' @param interval A length two vector containing the end-points of the
#' interval to be searched for the quantiles, passed to the uniroot function.
#' @param loc1,scale1,loc2,scale2 Location and scale parameters of the two
#' Gumbel distributions. The second location parameter must be greater than or
#' equal to the first location parameter.
#' @param log Logical; if \code{TRUE}, the log density is returned.
#' @param lower.tail Logical; if \code{TRUE} (default), probabilities are 
#' \verb{P[X <= x]}, otherwise, \verb{P[X > x]}
#' @param \dots Other arguments passed to uniroot.
#' @return \code{dgumbelx} gives the density function, \code{pgumbelx} gives
#' the distribution function, \code{qgumbelx} gives the quantile function, and
#' \code{rgumbelx} generates random deviates.
#' 
#' @seealso \code{\link[evd]{fgev}}, \code{\link{uniroot}}
#' 
#' @family topic.distributions
#' @concept continuous distribution
#' @concept extreme value theory
#' @concept extended Gumbel
#' @concept dpqr
#' 
#' @examples
#' 
#' dgumbelx(2:4, 0, 1.1, 1, 0.5)
#' pgumbelx(2:4, 0, 1.1, 1, 0.5)
#' qgumbelx(seq(0.9, 0.6, -0.1), interval = c(0,10), 0, 1.2, 2, 0.5)
#' rgumbelx(6, 0, 1.1, 1, 0.5)
#' p <- (1:9)/10
#' pgumbelx(qgumbelx(p, interval = c(0,10), 0, 0.5, 1, 2), 0, 0.5, 1, 2)
#' ## [1] 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9
#' 



#' @rdname dgumbelx
#' @export
"dgumbelx"<-
  function(x, loc1 = 0, scale1 = 1, loc2 = 0, scale2 = 1, log = FALSE)
  {
    if(min(scale1) < 0 || min(scale2) < 0) stop("invalid scale")
    if(any(loc1 > loc2)) stop("loc1 cannot be greater than loc2")
    x1 <- (x - loc1)/scale1
    x2 <- (x - loc2)/scale2
    d <- exp(-exp(-x1) + log(1/scale2) - x2 - exp(-x2)) + exp(-exp(-x2) + log(1/scale1) - x1 - exp(-x1)) 
    if(log) d <- log(d)
    d
  }


#' @rdname dgumbelx
#' @export
"pgumbelx"<-
  function(q, loc1 = 0, scale1 = 1, loc2 = 0, scale2 = 1, lower.tail = TRUE)
  {
    if(min(scale1) < 0 || min(scale2) < 0) stop("invalid scale")
    if(any(loc1 > loc2)) stop("loc1 cannot be greater than loc2")
    q1 <- (q - loc1)/scale1
    q2 <- (q - loc2)/scale2
    p <- exp(-exp(-q1)) * exp(-exp(-q2))
    if(!lower.tail) p <- 1 - p
    p
  }


#' @rdname dgumbelx
#' @export
"qgumbelx"<-
  function(p, interval, loc1 = 0, scale1 = 1, loc2 = 0, scale2 = 1, lower.tail = TRUE, ...)
  {
    if(min(p, na.rm = TRUE) <= 0 || max(p, na.rm = TRUE) >=1)
      stop("`p' must contain probabilities in (0,1)")
    if(min(scale1) < 0 || min(scale2) < 0) stop("invalid scale")
    if(any(loc1 > loc2)) stop("loc1 cannot be greater than loc2")
    if(!lower.tail) p <- 1 - p
    
    n <- length(p)
    out <- numeric(n)
    for(i in 1:n) {
      tmpfn <- function(z) exp(-(z - loc1)/scale1) + exp(-(z - loc2)/scale2) + log(p[i])
      out[i] <- uniroot(tmpfn, interval = interval, ...)$root
    }
    out
  }


#' @rdname dgumbelx
#' @export
"rgumbelx"<-
  function(n, loc1 = 0, scale1 = 1, loc2 = 0, scale2 = 1)
  {
    if(min(scale1) < 0 || min(scale2) < 0) stop("invalid scale")
    if(any(loc1 > loc2)) stop("loc1 cannot be greater than loc2")
    pmax(rgumbel(n = n, loc = loc1, scale = scale1), rgumbel(n = n, loc = loc2, scale = scale2))
  }
