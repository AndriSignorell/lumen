
#' Gompertz Distribution
#' 
#' The Gompertz distribution is a continuous distribution with a non-negative 
#' real support, commonly used to model human mortality and customer lifetime 
#' value. It is parameterized by a shape and a scale parameter, and is 
#' characterized by an exponentially increasing hazard rate.
#' 
#' The Gompertz distribution with `shape` parameter \eqn{a} and
#' `rate` parameter \eqn{b}{b} has probability density function
#' 
#' \deqn{f(x | a, b) = be^{ax}\exp(-b/a (e^{ax} - 1))}{f(x | a, b) = b exp(ax)
#' exp(-b/a (exp(ax) - 1))}
#' 
#' For \eqn{a=0} the Gompertz is equivalent to the exponential distribution
#' with constant hazard and rate \eqn{b}.
#' 
#' The probability distribution function is
#' \deqn{F(x | a, b) = 1 - \exp(-b/a (e^{ax} - 1))}{F(x | a, b) = 1 - exp(-b/a (exp(ax) - 1))}
#' 
#' Thus if \eqn{a} is negative, letting \eqn{x} tend to infinity shows that
#' there is a non-zero probability \eqn{1 - \exp(b/a)}{1 - exp(b/a)} of living
#' forever.  On these occasions `qgompertz()` and `rgompertz()` will
#' return `Inf`, and `pgompertz()` approaches
#' \eqn{1 - \exp(b/a)}{1 - exp(b/a)} rather than one.
#' 
#' A non-positive `rate` gives `NaN` with a warning, as in the base R
#' distribution functions.
#' 
#' **Note:** Some implementations of the Gompertz restrict \eqn{a} to be strictly
#' positive, which ensures that the probability of survival decreases to zero
#' as \eqn{x} increases to infinity.  The more flexible implementation given
#' here is consistent with `streg` in Stata.
#' 
#' The functions `dgompertz()` and similar available in the package
#' \pkg{eha} label the parameters the other way round, so that what is called
#' the `shape` there is called the `rate` here, and what is called
#' `1 / scale` there is called the `shape` here. The terminology here
#' is consistent with the exponential [dexp()] and Weibull
#' [dweibull()] distributions in R.
#' 
#' @name dpqr-gompertz
#' @aliases Gompertz dgompertz pgompertz qgompertz rgompertz
#' 
#' @param x,q vector of quantiles.
#' @param shape,rate vector of shape and rate parameters.
#' @param log,log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are \eqn{P(X
#' }{P(X <= x)}\eqn{\le x)}{P(X <= x)}, otherwise, \eqn{P(X > x)}{P(X > x)}.
#' @param p vector of probabilities.
#' @param n number of observations. If `length(n) > 1`, the length is
#' taken to be the number required.
#' @return `dgompertz()` gives the density, `pgompertz()` gives the
#' distribution function, `qgompertz()` gives the quantile function, and
#' `rgompertz()` generates random deviates.
#' 
#' @note
#' Based on code by Christopher Jackson previously published in
#' the \pkg{flexsurv} package, adapted to conform to package standards.
#'  
#' @seealso [distributions-overview]; [dexp()]
#' 
#' @references
#' Gompertz, B. (1825) On the nature of the function expressive of the law
#' of human mortality. *Philosophical Transactions of the Royal Society*,
#' **115**, 513--583.
#'
#' Stata Press (2007) *Stata Release 10 Manual: Survival Analysis
#' and Epidemiological Tables*. Stata Press.
#' 
#' @examples
#' 
#' dgompertz(1:3, shape = 0.1, rate = 0.2)
#' pgompertz(1:3, shape = 0.1, rate = 0.2)
#' qgompertz(seq(0.9, 0.6, -0.1), shape = 0.1, rate = 0.2)
#' rgompertz(6, shape = 0.1, rate = 0.2)
#' 
#' ## for shape = 0 the Gompertz reduces to the exponential distribution
#' all.equal(pgompertz(1:3, shape = 0, rate = 0.2), pexp(1:3, rate = 0.2))
#' 
#' ## a negative shape leaves a non-zero probability of living forever,
#' ## for which the quantile function returns Inf
#' qgompertz(0.9, shape = -0.5, rate = 0.2)
#' 
#' mgompertz(shape = 0.1, rate = 0.2)
#'  

#' @rdname dpqr-gompertz
#' @concept distribution-function
#' @concept demographics
#' @export
dgompertz <- function(x, shape, rate = 1, log = FALSE) {
  # the C++ side returns NaN for invalid parameters; the warning that goes
  # with it is raised here, once per call rather than once per element
  .checkGompertz(shape, rate)
  dgompertz_cpp(x, shape, rate, log)
}


#' @rdname dpqr-gompertz
#' @export
pgompertz <- function(q, shape, rate = 1, lower.tail = TRUE, log.p = FALSE) {
  .checkGompertz(shape, rate)
  pgompertz_cpp(q, shape, rate, lower.tail, log.p)
}


#' @rdname dpqr-gompertz
#' @export
qgompertz <- function(p, shape, rate = 1, lower.tail = TRUE, log.p = FALSE) {
  d     <- .dbase(.checkGompertz, lower.tail = lower.tail, log = log.p, 
                  p = p, shape = shape, rate = rate)
  ret   <- d$ret
  ind   <- d$ind
  p     <- d$p
  shape <- d$shape
  rate  <- d$rate
  
  s0 <- abs(shape) <= .gompertzShapeTol
  ret[ind][s0] <- qexp(p[s0], rate = rate[s0])
  sn0 <- !s0
  if (any(sn0)) {
    p     <- p[sn0]
    shape <- shape[sn0]
    rate  <- rate[sn0]
    asymp   <- 1 - exp(rate / shape)
    immortal <- shape < 0 & p > asymp
    ret[ind][sn0][immortal]  <- Inf
    ret[ind][sn0][!immortal] <- 1 / shape[!immortal] *
      log1p(-log1p(-p[!immortal]) * shape[!immortal] / rate[!immortal])
  }
  ret
}


#' @rdname dpqr-gompertz
#' @export
rgompertz <- function(n, shape, rate = 1) {
  r     <- .rbase(.checkGompertz, n = n, shape = shape, rate = rate)
  ret   <- r$ret
  ind   <- r$ind
  shape <- r$shape
  rate  <- r$rate
  
  ret[ind] <- qgompertz(p = runif(sum(ind)), shape = shape, rate = rate)
  ret
}


# == internal helper functions ===============================================


# parameter validity for the Gompertz distribution: any real shape is
# allowed, the rate must be positive (invalid entries yield NaN with a
# warning, matching the base R d/p/q/r convention). The test itself lives
# in C++ so that d, p, q and r cannot drift apart.
.checkGompertz <- function(shape, rate) {
  ok <- checkGompertz_cpp(shape, rate)
  if (any(!ok)) warning("Non-positive rate parameter")
  ok
}

# below this the shape counts as zero and the exponential limit is used;
# the same value is hard-wired in gompertz.cpp
.gompertzShapeTol <- 1e-12


### Standardised procedure for defining density, cumulative
### distribution, hazard and cumulative hazard functions for
### time-to-event distributions

.dbase <- function(checkFun, lower.tail=TRUE, log=FALSE, ...){
  args <- list(...)
  ## Vectorise all arguments, replicating to length of longest argument
  n <- max(sapply(args, length))
  for (i in seq_along(args)) {
    args[[i]] <- rep(args[[i]], length=n)
  }
  ret <- numeric(n)
  ## Check for parameters out of range, give warning and return NaN
  ## for those
  check.ret <- do.call(checkFun, args[-1])
  ret[!check.ret] <- NaN
  for (i in seq_along(args))
    ret[is.nan(args[[i]])] <- NaN
  ## name of first arg is x for PDF, haz, or cum haz, q for CDF and p for quantile function
  stopifnot( !(names(args)[1]=="x" && lower.tail==FALSE))
  if (names(args)[1] %in% c("x","q")){
    x <- args[[1]]
    ## PDF, CDF, hazard and cumulative hazard is 0 for any negative time
    ret[!is.nan(ret) & (x<0)] <- if (lower.tail) { if (log) -Inf else 0 } else { if (log) 0 else 1 }
  }
  if (names(args)[1] == "p") {
    p <- args[[1]]
    if (log) p <- exp(p)
    if (!lower.tail) p <- 1 - p
    args[[1]] <- p
    ret[p < 0 | p > 1] <- NaN
    ## should be 0,Inf for p=0,1, but hopefully always handled anyway
    ## Result is NA if x or a parameter is NA
  }
  ## Result is NA if x or a parameter is NA
  nas <- rep(FALSE, n)
  for (i in seq_along(args)) nas <- nas | (is.na(args[[i]]) & !is.nan(args[[i]]))
  ret[nas] <- NA
  ind <- !is.nan(ret) & !nas
  if (names(args)[1] %in% c("x", "q")) ind <- ind & (x>=0)
  ## Any remaining elements of vector are filled in by standard
  ## formula for hazard
  li <- list(ret=ret, ind=ind)
  for(i in seq_along(args)) args[[i]] <- args[[i]][ind]
  c(li, args)
}

### Standardised procedure for defining random sampling functions

.rbase <- function(checkFun, n, ...){
  ## Vectorise all arguments, replicating to sample length
  if (length(n) > 1) n <- length(n)
  args <- list(...)
  for (i in seq_along(args)) {
    args[[i]] <- rep(args[[i]], length=n)
  }
  ret <- numeric(n)
  ## Check for parameters out of range, give warning and return NaN
  ## for those
  check.ret <- do.call(checkFun, args)
  ret[!check.ret] <- NaN
  for (i in seq_along(args))
    ret[is.nan(args[[i]])] <- NaN
  nas <- rep(FALSE, n)
  for (i in seq_along(args)) nas <- nas | (is.na(args[[i]]) & !is.nan(args[[i]]))
  ret[nas] <- NA
  ind <- !is.nan(ret) & !nas
  li <- list(ret=ret, ind=ind)
  for(i in seq_along(args)) args[[i]] <- args[[i]][ind]
  c(li, args)
}
