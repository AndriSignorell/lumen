

#' Distributions of Order Statistics
#' 
#' Density, distribution, and random generation functions for a selected 
#' order statistic (the j-th largest or smallest value) from a sample of 
#' a given size drawn from any specified distribution, derived analytically 
#' using the beta distribution representation of order statistics.
#' 
#' @name dpqr-order
#' @aliases dorder porder rorder
#' 
#' @param x,q vector of quantiles.
#' @param n number of observations.
#' @param dFun,pFun,qFun density, distribution and quantile function
#' of the specified distribution. The density function must have a `log`
#' argument (a simple wrapper can always be constructed to achieve this).
#' @param \dots parameters of the specified distribution.
#' @param distn a character string, optionally specified as an alternative to
#' `dFun`, `pFun` and `qFun` such that the density,
#' distribution and quantile functions are formed upon the addition of the
#' prefixes `d`, `p` and `q` respectively.
#' @param mlen the number of independent variables.
#' @param j the order statistic, taken as the `j`th largest (default) or
#' smallest of `mlen`, according to the value of `largest`.
#' @param largest logical; if `TRUE` (default) use the `j`th largest
#' order statistic, otherwise use the `j`th smallest.
#' @param log,log.p logical; if `TRUE`, probabilities `p` are given as
#' `log(p)` and the density is returned on the log scale.
#' @param lower.tail logical; if `TRUE` (default) probabilities are 
#' \verb{P[X <= x]}, otherwise P\verb{[X > x]}.
#' @return `dorder()` gives the density function and `porder()`
#' gives the distribution function of a selected order statistic from a
#' sample of size `mlen`, from a specified distribution.
#' `rorder()` generates random deviates. There is no quantile function
#' for order statistics (`qorder()` does not exist).
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


#' @rdname dpqr-order
#' @export
dorder <- function(x, dFun, pFun, ..., distn, mlen = 1, j = 1, largest = TRUE,
           log = FALSE)
  {
    .checkOrderIndex(mlen, j)
    if(!largest) j <- mlen + 1 - j
    if(missing(dFun))
      dFun <- get(paste0("d", distn), mode="function")
    if(missing(pFun))
      pFun <- get(paste0("p", distn), mode="function")
    dens <- dFun(x, ..., log = TRUE)
    ok <- !is.infinite(dens)
    Fx <- pFun(x, ...)[ok]
    # each exponent vanishes at one end of the support, where the
    # corresponding logarithm is -Inf
    lFx <- (if(mlen == j) 0 else (mlen-j) * log(Fx)) +
           (if(j == 1L)   0 else (j-1) * log1p(-Fx))
    comb <- lgamma(mlen+1) - lgamma(j) - lgamma(mlen-j+1)
    d <- numeric(length(x))
    d[ok]  <- comb + dens[ok] + lFx
    d[!ok] <- -Inf
    if(!log) d <- exp(d)
    d
  }




#' @rdname dpqr-order
#' @export
porder <- function(q, pFun, ..., distn, mlen = 1, j = 1, largest = TRUE,
           lower.tail = TRUE, log.p = FALSE)
  {
    .checkOrderIndex(mlen, j)
    if(largest) svec <- (mlen+1-j):mlen
    else  svec <- 0:(j-1)
    if(missing(pFun))
      pFun <- get(paste0("p", distn), mode="function")
    Fx <- pFun(q, ...)
    store <- matrix(0, nrow = length(q), ncol = j)
    for(k in 1:j) {
      s <- svec[k]
      # a zero exponent cancels the -Inf of log(0) at the ends of the support
      store[,k] <- exp(lchoose(mlen, s) +
                       (if(s == 0)    0 else s * log(Fx)) +
                       (if(s == mlen) 0 else (mlen-s) * log1p(-Fx)))
    }
    p <- rowSums(store)
    if(largest != lower.tail) p <- 1 - p
    if(log.p) log(p) else p
  }


#' @rdname dpqr-order
#' @export
rorder <- function(n, qFun, ..., distn,  mlen = 1, j = 1, largest = TRUE)
  {
    .checkOrderIndex(mlen, j)
    if(!largest) j <- mlen+1-j
    if(missing(qFun))
      qFun <- get(paste0("q", distn), mode="function")
    qFun(rbeta(n, mlen+1-j, j), ...)
  }




# == internal helper functions ===============================================

#' @noRd
.checkOrderIndex <- function(mlen, j = 1) {
  .assertScalar(mlen, lower = 1, integerValued = TRUE)
  .assertScalar(j,    lower = 1, integerValued = TRUE)
  if(j > mlen)
    stop("'j' cannot be greater than 'mlen'", call. = FALSE)
}
