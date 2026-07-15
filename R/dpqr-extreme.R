
#' Distributions of Maxima and Minima
#' 
#' Density, distribution, quantile, and random generation functions for 
#' the maximum or minimum of a given number of independent and identically 
#' distributed random variables from any specified distribution, derived 
#' analytically from the underlying distribution function.
#' 
#' 
#' Density function, distribution function, quantile function and random
#' generation for the maximum/minimum of a given number of independent
#' variables from a specified distribution.
#' 
#' 
#' @name dpqr-extreme
#' @aliases dextreme pextreme qextreme rextreme
#' 
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param n number of observations
#' @param densfun,distnfun,quantfun density, distribution and quantile function
#' of the specified distribution. The density function must have a \code{log}
#' argument (a simple wrapper can always be constructed to achieve this).
#' @param \dots parameters of the specified distribution
#' @param distn a character string, optionally given as an alternative to
#' \code{densfun}, \code{distnfun} and \code{quantfun} such that the density,
#' distribution and quantile functions are formed upon the addition of the
#' prefixes \code{d}, \code{p} and \code{q} respectively
#' @param mlen the number of independent variables
#' @param largest logical; if \code{TRUE} (default) use maxima, otherwise
#' minima
#' @param log logical; if \code{TRUE}, the log density is returned
#' @param lower.tail logical; if \code{TRUE} (default) probabilities are 
#' \verb{P[X <= x]}, otherwise P\verb{[X > x]}
#' 
#' @return \code{dextreme()} gives the density function, \code{pextreme()}
#' gives the distribution function and \code{qextreme()} gives the quantile
#' function of the maximum/minimum of \code{mlen} independent variables from a
#' specified distribution. \code{rextreme()} generates random deviates.
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
#' dextreme(2:4, dnorm, pnorm, mean = 0.5, sd = 1.2, mlen = 5)
#' dextreme(2:4, distn = "norm", mean = 0.5, sd = 1.2, mlen = 5)
#' dextreme(2:4, distn = "exp", mlen = 2, largest = FALSE)
#' pextreme(2:4, distn = "exp", rate = 1.2, mlen = 2)
#' qextreme(seq(0.9, 0.6, -0.1), distn = "exp", rate = 1.2, mlen = 2)
#' rextreme(5, qgamma, shape = 1, mlen = 10)
#' p <- (1:9)/10
#' pexp(qextreme(p, distn = "exp", rate = 1.2, mlen = 1), rate = 1.2)
#' ## [1] 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9
#' 
#' 
#' @rdname dpqr-extreme
#' @export
dextreme <- function(x, densfun, distnfun, ..., distn, mlen = 1, largest = TRUE, log = FALSE)
  {
    .checkOrderIndex(mlen)
    if(missing(densfun))
      densfun <- get(paste("d", distn, sep=""), mode="function")
    if(missing(distnfun))
      distnfun <- get(paste("p", distn, sep=""), mode="function")
    dens <- densfun(x, ..., log = TRUE)
    distn <- distnfun(x, ...)[!is.infinite(dens)]
    if(!largest) distn <- 1 - distn
    distn <- (mlen-1) * log(distn)
    d <- numeric(length(x))
    d[!is.infinite(dens)] <- log(mlen) + dens[!is.infinite(dens)] + distn
    d[is.infinite(dens)] <- -Inf
    if(!log) d <- exp(d)
    d
  }


#' @rdname dpqr-extreme
#' @export
pextreme <- function(q, distnfun, ..., distn, mlen = 1, largest = TRUE, lower.tail = TRUE)
  {
    .checkOrderIndex(mlen)
    if(missing(distnfun))
      distnfun <- get(paste("p", distn, sep=""), mode="function")
    distn <- distnfun(q, ...)
    if(!largest) distn <- 1-distn
    p <- distn^mlen
    if(largest != lower.tail) p <- 1 - p
    p
  }


#' @rdname dpqr-extreme
#' @export
qextreme <- function(p, quantfun, ..., distn, mlen = 1, largest = TRUE, lower.tail = TRUE)
  {
    if(min(p, na.rm = TRUE) <= 0 || max(p, na.rm = TRUE) >=1)
      stop("`p' must contain probabilities in (0,1)")
    .checkOrderIndex(mlen)
    if(missing(quantfun))
      quantfun <- get(paste("q", distn, sep=""), mode="function")
    if(!lower.tail) p <- 1 - p
    if(largest) 
      quantfun(p^(1/mlen), ...)
    else
      quantfun(1-(1-p)^(1/mlen), ...)
  }


#' @rdname dpqr-extreme
#' @export
rextreme <- function(n, quantfun, ..., distn, mlen = 1, largest = TRUE)
  {
    .checkOrderIndex(mlen)
    if(missing(quantfun))
      quantfun <- get(paste("q", distn, sep=""), mode="function")
    if(largest)
      quantfun(rbeta(n, mlen, 1), ...)
    else
      quantfun(rbeta(n, 1, mlen), ...) 
  }

