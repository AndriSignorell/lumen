
#' Distributions of Order Statistics
#' 
#' Density, distribution, and random generation functions for a selected 
#' order statistic (the j-th largest or smallest value) from a sample of 
#' a given size drawn from any specified distribution, derived analytically 
#' using the beta distribution representation of order statistics.
#' 
#' Density function, distribution function and random generation for a selected
#' order statistic of a given number of independent variables from a specified
#' distribution.
#' 
#' 
#' @name dpqr-order
#' @aliases dorder porder rorder
#' @param x,q vector of quantiles
#' @param n number of observations
#' @param dFun,pFun,qFun density, distribution and quantile function
#' of the specified distribution. The density function must have a \code{log}
#' argument (a simple wrapper can always be constructed to achieve this).
#' @param \dots parameters of the specified distribution
#' @param distn a character string, optionally specified as an alternative to
#' \code{dFun}, \code{pFun} and \code{qFun} such that the density,
#' distribution and quantile functions are formed upon the addition of the
#' prefixes \code{d}, \code{p} and \code{q} respectively
#' @param mlen the number of independent variables
#' @param j the order statistic, taken as the \code{j}th largest (default) or
#' smallest of \code{mlen}, according to the value of \code{largest}
#' @param largest logical; if \code{TRUE} (default) use the \code{j}th largest
#' order statistic, otherwise use the \code{j}th smallest
#' @param log logical; if \code{TRUE}, the log density is returned
#' @param lower.tail logical; if \code{TRUE} (default) probabilities are 
#' \verb{P[X <= x]}, otherwise P\verb{[X > x]}
#' @return \code{dorder()} gives the density function and \code{porder()}
#' gives the distribution function of a selected order statistic from a
#' sample of size \code{mlen}, from a specified distribution.
#' \code{rorder()} generates random deviates. There is no quantile function
#' for order statistics (\code{qorder()} does not exist).
#' 
#' @note
#' Based on code by Alec Stephenson previously published in
#' the \pkg{evd} package, adapted to conform to package standards.
#' 
#' @seealso [distributions-overview]
#' @concept distribution-function
#' @concept order-statistic
#' 
#' @examples
#' 
#' dorder(2:4, dnorm, pnorm, mean = 0.5, sd = 1.2, mlen = 5, j = 2)
#' dorder(2:4, distn = "norm", mean = 0.5, sd = 1.2, mlen = 5, j = 2)
#' dorder(2:4, distn = "exp", mlen = 2, j = 2)
#' porder(2:4, distn = "exp", rate = 1.2, mlen = 2, j = 2)
#' rorder(5, qgamma, shape = 1, mlen = 10, j = 2)
#' 

#' @noRd
.checkOrderIndex <- function(mlen, j = 1) {
  if(!is.numeric(mlen) || length(mlen) != 1 || mlen < 1 ||
     mlen %% 1 != 0) 
    stop("argument 'mlen' must be a positive integer")
  if(!is.numeric(j) || length(j) != 1 || j < 1 || j %% 1 != 0) 
    stop("argument 'j' must be a positive integer")
  if(j > mlen)
    stop("argument 'j' cannot be greater than 'mlen'")
}


#' @rdname dpqr-order
#' @export
dorder <- function(x, dFun, pFun, ..., distn, mlen = 1, j = 1, largest = TRUE,
           log = FALSE)
  {
    .checkOrderIndex(mlen, j)
    if(!largest) j <- mlen + 1 - j
    if(missing(dFun))
      dFun <- get(paste("d", distn, sep=""), mode="function")
    if(missing(pFun))
      pFun <- get(paste("p", distn, sep=""), mode="function")
    dens <- dFun(x, ..., log = TRUE)
    Fx <- pFun(x, ...)[!is.infinite(dens)]
    Fx <- (mlen-j) * log(Fx) + (j-1) * log(1-Fx)
    comb <- lgamma(mlen+1) - lgamma(j) - lgamma(mlen-j+1)
    d <- numeric(length(x))
    d[!is.infinite(dens)] <- comb + dens[!is.infinite(dens)] + Fx
    d[is.infinite(dens)] <- -Inf
    if(!log) d <- exp(d)
    d
  }




#' @rdname dpqr-order
#' @export
porder <- function(q, pFun, ..., distn, mlen = 1, j = 1, largest = TRUE,
           lower.tail = TRUE)
  {
    .checkOrderIndex(mlen, j)
    if(largest) svec <- (mlen+1-j):mlen
    else  svec <- 0:(j-1)
    if(missing(pFun))
      pFun <- get(paste("p", distn, sep=""), mode="function")
    Fx <- pFun(q, ...)
    store <- matrix(0,nrow=length(q),ncol=j)
    for(k in 1:j)
      store[,k] <- exp(lchoose(mlen, svec[k]) + svec[k]*log(Fx) +
                         (mlen-svec[k])*log(1-Fx))
    p <- apply(store,1,sum)
    if(largest != lower.tail) p <- 1 - p
    p
  }


#' @rdname dpqr-order
#' @export
rorder <- function(n, qFun, ..., distn,  mlen = 1, j = 1, largest = TRUE)
  {
    .checkOrderIndex(mlen, j)
    if(!largest) j <- mlen+1-j
    if(missing(qFun))
      qFun <- get(paste("q", distn, sep=""), mode="function")
    qFun(rbeta(n, mlen+1-j, j), ...)
  }

