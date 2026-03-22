
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
#' @name dorder
#' @aliases dorder porder rorder
#' @param x,q Vector of quantiles.
#' @param n Number of observations.
#' @param densfun,distnfun,quantfun Density, distribution and quantile function
#' of the specified distribution. The density function must have a \code{log}
#' argument (a simple wrapper can always be constructed to achieve this).
#' @param \dots Parameters of the specified distribution.
#' @param distn A character string, optionally specified as an alternative to
#' \code{densfun}, \code{distnfun} and \code{quantfun} such that the density,
#' distribution and quantile functions are formed upon the addition of the
#' prefixes \code{d}, \code{p} and \code{q} respectively.
#' @param mlen The number of independent variables.
#' @param j The order statistic, taken as the \code{j}th largest (default) or
#' smallest of \code{mlen}, according to the value of \code{largest}.
#' @param largest Logical; if \code{TRUE} (default) use the \code{j}th largest
#' order statistic, otherwise use the \code{j}th smallest.
#' @param log Logical; if \code{TRUE}, the log density is returned.
#' @param lower.tail Logical; if \code{TRUE} (default) probabilities are 
#' \verb{P[X <= x]}, otherwise P\verb{[X > x]}.
#' @return \code{dorder} gives the density function, \code{porder} gives the
#' distribution function and \code{qorder} gives the quantile function of a
#' selected order statistic from a sample of size \code{mlen}, from a specified
#' distribution. \code{rorder} generates random deviates.
#' 
#' @note
#' Based on code by Alec Stephenson. 
#' 
#' @family topic.extremevalue.distributions
#' @concept continuous distribution
#' @concept order statistics
#' @concept ranks
#' 
#' @examples
#' 
#' dorder(2:4, dnorm, pnorm, mean = 0.5, sd = 1.2, mlen = 5, j = 2)
#' dorder(2:4, distn = "norm", mean = 0.5, sd = 1.2, mlen = 5, j = 2)
#' dorder(2:4, distn = "exp", mlen = 2, j = 2)
#' porder(2:4, distn = "exp", rate = 1.2, mlen = 2, j = 2)
#' rorder(5, qgamma, shape = 1, mlen = 10, j = 2)
#' 

#' @rdname dorder
#' @export
"dorder"<-
  function(x, densfun, distnfun, ..., distn, mlen = 1, j = 1, largest = TRUE,
           log = FALSE)
  {
    if(!is.numeric(mlen) || length(mlen) != 1 || mlen < 1 ||
       mlen %% 1 != 0) 
      stop("`mlen' must be a non-negative integer")
    if(!is.numeric(j) || length(j) != 1 || j < 1 || j %% 1 != 0) 
      stop("`j' must be a non-negative integer")
    if(j > mlen)
      stop("`j' cannot be greater than `mlen'")
    if(!largest) j <- mlen + 1 - j
    if(missing(densfun))
      densfun <- get(paste("d", distn, sep=""), mode="function")
    if(missing(distnfun))
      distnfun <- get(paste("p", distn, sep=""), mode="function")
    dens <- densfun(x, ..., log = TRUE)
    distn <- distnfun(x, ...)[!is.infinite(dens)]
    distn <- (mlen-j) * log(distn) + (j-1) * log(1-distn)
    comb <- lgamma(mlen+1) - lgamma(j) - lgamma(mlen-j+1)
    d <- numeric(length(x))
    d[!is.infinite(dens)] <- comb + dens[!is.infinite(dens)] + distn
    d[is.infinite(dens)] <- -Inf
    if(!log) d <- exp(d)
    d
  }




#' @rdname dorder
#' @export
"porder"<-
  function(q, distnfun, ..., distn, mlen = 1, j = 1, largest = TRUE,
           lower.tail = TRUE)
  {
    if(!is.numeric(mlen) || length(mlen) != 1 || mlen < 1 ||
       mlen %% 1 != 0) 
      stop("`mlen' must be a non-negative integer")
    if(!is.numeric(j) || length(j) != 1 || j < 1 || j %% 1 != 0) 
      stop("`j' must be a non-negative integer")
    if(j > mlen)
      stop("`j' cannot be greater than `mlen'")
    lachooseb <- function(a,b) lgamma(a+1) - lgamma(b+1) - lgamma(a-b+1)
    if(largest) svec <- (mlen+1-j):mlen
    else  svec <- 0:(j-1)
    if(missing(distnfun))
      distnfun <- get(paste("p", distn, sep=""), mode="function")
    distn <- distnfun(q, ...)
    store <- matrix(0,nrow=length(q),ncol=j)
    for(k in 1:j)
      store[,k] <- exp(lachooseb(mlen,svec[k]) + svec[k]*log(distn) +
                         (mlen-svec[k])*log(1-distn))
    p <- apply(store,1,sum)
    if(largest != lower.tail) p <- 1 - p
    p
  }


#' @rdname dorder
#' @export
"rorder"<-
  function(n, quantfun, ..., distn,  mlen = 1, j = 1, largest = TRUE)
  {
    if(!is.numeric(mlen) || length(mlen) != 1 || mlen < 1 ||
       mlen %% 1 != 0) 
      stop("`mlen' must be a non-negative integer")
    if(!is.numeric(j) || length(j) != 1 || j < 1 || j %% 1 != 0) 
      stop("`j' must be a non-negative integer")
    if(j > mlen)
      stop("`j' cannot be greater than `mlen'")
    if(!largest) j <- mlen+1-j
    if(missing(quantfun))
      quantfun <- get(paste("q", distn, sep=""), mode="function")
    quantfun(rbeta(n, mlen+1-j, j), ...)
  }

