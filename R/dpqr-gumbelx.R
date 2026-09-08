
#' Maxima of Two Gumbel Distributions
#' 
#' The extended Gumbel distribution models the maximum of two independent 
#' Gumbel-distributed random variables with potentially different location 
#' and scale parameters. It is parameterized by two pairs of location and 
#' scale parameters.
#' 
#' Density function, distribution function, quantile function and random
#' generation for the maxima of two Gumbel distributions, each with different
#' location and scale parameters.
#' 
#' 
#' @name dpqr-gumbelx
#' @aliases dgumbelx pgumbelx qgumbelx rgumbelx
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations.
#' @param interval a length two vector containing the end-points of the
#' interval to be searched for the quantiles, passed to [uniroot()]. By
#' default a bracketing interval is derived from the quantiles of the two
#' Gumbel margins.
#' @param loc1,scale1,loc2,scale2 location and scale parameters of the two
#' Gumbel distributions. The distribution is symmetric in the two margins,
#' so their order is immaterial.
#' @param log,log.p logical; if `TRUE`, probabilities `p` are given as
#' `log(p)` and the density is returned on the log scale.
#' @param lower.tail logical; if `TRUE` (default), probabilities are 
#' \verb{P[X <= x]}, otherwise, \verb{P[X > x]}.
#' @param \dots other arguments passed to uniroot.
#' @return `dgumbelx()` gives the density function, `pgumbelx()`
#' gives the distribution function, `qgumbelx()` gives the quantile
#' function, and `rgumbelx()` generates random deviates.
#' 
#' @note
#' Based on code by Alec Stephenson previously published in
#' the \pkg{evd} package, adapted to conform to package standards.
#' 
#' @seealso [distributions-overview]; [uniroot()], which
#' `qgumbelx()` uses for root finding
#' 
#' @concept distribution-function
#' @concept extreme-value
#' 
#' @examples
#' 
#' dgumbelx(2:4, 0, 1.1, 1, 0.5)
#' pgumbelx(2:4, 0, 1.1, 1, 0.5)
#' qgumbelx(seq(0.9, 0.6, -0.1), 0, 1.2, 2, 0.5)
#' rgumbelx(6, 0, 1.1, 1, 0.5)
#' p <- (1:9)/10
#' pgumbelx(qgumbelx(p, 0, 0.5, 1, 2), 0, 0.5, 1, 2)
#' ## [1] 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9
#' 



#' @rdname dpqr-gumbelx
#' @export
dgumbelx <- function(x, loc1 = 0, scale1 = 1, loc2 = 0, scale2 = 1,
                     log = FALSE)
  {
    .assertPositive(scale1)
    .assertPositive(scale2)
    x1 <- (x - loc1)/scale1
    x2 <- (x - loc2)/scale2

    # f1(x) F2(x) + f2(x) F1(x), summed on the log scale
    l1 <- -exp(-x1) - log(scale2) - x2 - exp(-x2)
    l2 <- -exp(-x2) - log(scale1) - x1 - exp(-x1)
    hi <- pmax(l1, l2)
    d  <- hi + log1p(exp(-abs(l1 - l2)))
    d[is.infinite(hi) & hi < 0] <- -Inf

    if(log) d else exp(d)
  }


#' @rdname dpqr-gumbelx
#' @export
pgumbelx <- function(q, loc1 = 0, scale1 = 1, loc2 = 0, scale2 = 1,
                     lower.tail = TRUE, log.p = FALSE)
  {
    .assertPositive(scale1)
    .assertPositive(scale2)
    q1 <- (q - loc1)/scale1
    q2 <- (q - loc2)/scale2

    lp <- -exp(-q1) - exp(-q2)
    if(lower.tail) {
      if(log.p) lp else exp(lp)
    } else {
      if(log.p) log(-expm1(lp)) else -expm1(lp)
    }
  }


#' @rdname dpqr-gumbelx
#' @export
qgumbelx <- function(p, loc1 = 0, scale1 = 1, loc2 = 0, scale2 = 1,
                     lower.tail = TRUE, log.p = FALSE, interval = NULL, ...)
  {
    .assertScalar(loc1);   .assertScalar(scale1, lower = 0, strictLower = TRUE)
    .assertScalar(loc2);   .assertScalar(scale2, lower = 0, strictLower = TRUE)
    p <- .qProb(p, lower.tail = lower.tail, log.p = log.p)

    # uniroot's default tolerance would leave the quantiles accurate to
    # about 1e-4 only
    dots <- list(...)
    if(is.null(dots$tol)) dots$tol <- .Machine$double.eps^0.5

    out <- numeric(length(p))
    for(i in seq_along(p)) {

      if(is.na(p[i])) { out[i] <- p[i]; next }
      if(p[i] == 0)   { out[i] <- -Inf;  next }
      if(p[i] == 1)   { out[i] <-  Inf;  next }

      # F = F1 F2 is bounded by min(F1, F2) from above and, for z beyond
      # both sqrt(p) margin quantiles, by p from below -- an exact bracket
      lo <- max(qgumbel(p[i], loc1, scale1), qgumbel(p[i], loc2, scale2))
      hi <- max(qgumbel(sqrt(p[i]), loc1, scale1),
                qgumbel(sqrt(p[i]), loc2, scale2))

      tmpfn <- function(z)
        exp(-(z - loc1)/scale1) + exp(-(z - loc2)/scale2) + log(p[i])

      out[i] <- do.call(uniroot,
                        c(list(f = tmpfn,
                               interval = if(is.null(interval)) c(lo, hi)
                                          else interval), dots))$root
    }
    out
  }


#' @rdname dpqr-gumbelx
#' @export
rgumbelx <- function(n, loc1 = 0, scale1 = 1, loc2 = 0, scale2 = 1)
  {
    .assertPositive(scale1)
    .assertPositive(scale2)
    pmax(rgumbel(n = n, loc = loc1, scale = scale1), rgumbel(n = n, loc = loc2, scale = scale2))
  }
