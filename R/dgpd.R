

#' The Generalized Pareto Distribution
#' 
#' Density function, distribution function, quantile function and random
#' generation for the generalized Pareto distribution (GPD) with location,
#' scale and shape parameters.
#' 
#' The generalized Pareto distribution function (Pickands, 1975) with
#' parameters \eqn{\code{loc} = a}, \eqn{\code{scale} = b} and
#' \eqn{\code{shape} = s} is \deqn{G(z) = 1 - \{1+s(z-a)/b\}^{-1/s}}{ G(z) = 1
#' - {1+s(z-a)/b}^(-1/s)} for \eqn{1+s(z-a)/b > 0} and \eqn{z > a}, where
#' \eqn{b > 0}.  If \eqn{s = 0} the distribution is defined by continuity.
#' 
#' @name dgpd
#' @aliases dgpd pgpd qgpd rgpd
#' @param x,q Vector of quantiles.
#' @param p Vector of probabilities.
#' @param n Number of observations.
#' @param loc,scale,shape Location, scale and shape parameters; the
#' \code{shape} argument cannot be a vector (must have length one).
#' @param log Logical; if \code{TRUE}, the log density is returned.
#' @param lower.tail Logical; if \code{TRUE} (default), probabilities are 
#' \verb{P[X <= x]}, otherwise, P\verb{[X > x]}
#' @return \code{dgpd} gives the density function, \code{pgpd} gives the
#' distribution function, \code{qgpd} gives the quantile function, and
#' \code{rgpd} generates random deviates.
#' 
#' @seealso \code{\link[evd]{fpot}}
#' 
#' @references Pickands, J. (1975) Statistical inference using extreme order
#' statistics.  \emph{Annals of Statistics}, \bold{3}, 119--131.
#' 
#' @family topic.distributions
#' @concept continuous distribution
#' @concept extreme value theory
#' @concept peaks over threshold
#' @concept GPD
#' @concept dpqr
#' 
#' @examples
#' 
#' dgpd(2:4, 1, 0.5, 0.8)
#' pgpd(2:4, 1, 0.5, 0.8)
#' qgpd(seq(0.9, 0.6, -0.1), 2, 0.5, 0.8)
#' rgpd(6, 1, 0.5, 0.8)
#' p <- (1:9)/10
#' pgpd(qgpd(p, 1, 2, 0.8), 1, 2, 0.8)
#' ## [1] 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9
#' 

#' @rdname dgpd
#' @export
"dgpd"<-
  function(x, loc = 0, scale = 1, shape = 0, log = FALSE)
  {
    if(min(scale) <= 0) stop("invalid scale")
    if(length(shape) != 1) stop("invalid shape")
    d <- (x - loc)/scale
    nn <- length(d)
    scale <- rep(scale, length.out = nn)
    index <- (d > 0 & ((1 + shape * d) > 0)) | is.na(d)
    if(shape == 0) {
      d[index] <- log(1/scale[index]) - d[index]
      d[!index] <- -Inf
    }
    else {
      d[index] <- log(1/scale[index]) - (1/shape + 1) *
        log(1 + shape * d[index])
      d[!index] <- -Inf
    }
    if(!log) d <- exp(d)
    d
  }


#' @rdname dgpd
#' @export
"pgpd" <-
  function(q, loc = 0, scale = 1, shape = 0, lower.tail = TRUE)
  {
    if(min(scale) <= 0) stop("invalid scale")
    if(length(shape) != 1) stop("invalid shape")
    q <- pmax(q - loc, 0)/scale
    if(shape == 0) p <- 1 - exp(-q)
    else {
      p <- pmax(1 + shape * q, 0)
      p <- 1 - p^(-1/shape)
    }
    if(!lower.tail) p <- 1 - p
    p
  }

#' @rdname dgpd
#' @export
"qgpd"<-
  function(p, loc = 0, scale = 1, shape = 0, lower.tail = TRUE)
  {
    if(min(p, na.rm = TRUE) <= 0 || max(p, na.rm = TRUE) >=1)
      stop("`p' must contain probabilities in (0,1)")
    if(min(scale) < 0) stop("invalid scale")
    if(length(shape) != 1) stop("invalid shape")
    if(lower.tail) p <- 1 - p
    if(shape == 0) return(loc - scale*log(p))
    else return(loc + scale * (p^(-shape) - 1) / shape)
  }


#' @rdname dgpd
#' @export
"rgpd"<-
  function(n, loc = 0, scale = 1, shape = 0)
  {
    if(min(scale) < 0) stop("invalid scale")
    if(length(shape) != 1) stop("invalid shape")
    if(shape == 0) return(loc + scale*rexp(n))
    else return(loc + scale * (runif(n)^(-shape) - 1) / shape)
  }


