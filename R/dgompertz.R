
#' The Gompertz distribution
#' 
#' The Gompertz distribution is a continuous distribution with a non-negative 
#' real support, commonly used to model human mortality and customer lifetime 
#' value. It is parameterized by a shape and a scale parameter, and is 
#' characterized by an exponentially increasing hazard rate.
#' 
#' The Gompertz distribution with \code{shape} parameter \eqn{a} and
#' \code{rate} parameter \eqn{b}{b} has probability density function
#' 
#' \deqn{f(x | a, b) = be^{ax}\exp(-b/a (e^{ax} - 1))}{f(x | a, b) = b exp(ax)
#' exp(-b/a (exp(ax) - 1))}
#' 
#' For \eqn{a=0} the Gompertz is equivalent to the exponential distribution
#' with constant hazard and rate \eqn{b}.
#' 
#' The probability distribution function is \deqn{F(x | a, b) = 1 - \exp(-b/a
#' }{F(x | a, b) = 1 - exp(-b/a (exp(ax) - 1))}\deqn{(e^{ax} - 1))}{F(x | a, b)
#' = 1 - exp(-b/a (exp(ax) - 1))}
#' 
#' Thus if \eqn{a} is negative, letting \eqn{x} tend to infinity shows that
#' there is a non-zero probability \eqn{1 - \exp(b/a)}{1 - exp(b/a)} of living
#' forever.  On these occasions \code{qgompertz} and \code{rgompertz} will
#' return \code{Inf}.
#' 
#' \strong{Note:} \verb{    } Some implementations of the Gompertz restrict \eqn{a} to be strictly
#' positive, which ensures that the probability of survival decreases to zero
#' as \eqn{x} increases to infinity.  The more flexible implementation given
#' here is consistent with \code{streg} in Stata.
#' 
#' The functions \code{dgompertz} and similar available in the package
#' \pkg{eha} label the parameters the other way round, so that what is called
#' the \code{shape} there is called the \code{rate} here, and what is called
#' \code{1 / scale} there is called the \code{shape} here. The terminology here
#' is consistent with the exponential \code{\link{dexp}} and Weibull
#' \code{\link{dweibull}} distributions in R.
#' 
#' @name dgompertz
#' @aliases Gompertz dgompertz pgompertz qgompertz rgompertz
#' 
#' @param x,q vector of quantiles.
#' @param shape,rate vector of shape and rate parameters.
#' @param log,log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are \eqn{P(X
#' }{P(X <= x)}\eqn{\le x)}{P(X <= x)}, otherwise, \eqn{P(X > x)}{P(X > x)}.
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is
#' taken to be the number required.
#' @return \code{dgompertz} gives the density, \code{pgompertz} gives the
#' distribution function, \code{qgompertz} gives the quantile function, and
#' \code{rgompertz} generates random deviates.
#' 
#' @note
#' Adapted from \code{flexsurv::qgompertz()} by Christopher Jackson
#' to conform to package standards.
#'  
#' @seealso \code{\link{dexp}}
#' 
#' @references
#' Gompertz, B. (1825) On the nature of the function expressive of the law
#' of human mortality. \emph{Philosophical Transactions of the Royal Society},
#' \bold{115}, 513--583.
#'
#' Stata Press (2007) \emph{Stata Release 10 Manual: Survival Analysis
#' and Epidemiological Tables}. Stata Press.
#'  

#' @rdname dgompertz
#' @family dist.other
#' @concept distributions
#' @concept survival-analysis
#'
#'
#' @export
dgompertz <- function(x, shape, rate=1, log=FALSE) {
  # this is a verbatim copy from the package flexsurv (Christopher Jackson)
  dgompertz_cpp(x, shape, rate, log)
}


#' @rdname dgompertz
#' @export
pgompertz <- function(q, shape, rate=1, lower.tail = TRUE, log.p = FALSE) {
  pgompertz_cpp(q, shape, rate, lower.tail, log.p)
}


#' @rdname dgompertz
#' @export
qgompertz <- function(p, shape, rate = 1, lower.tail = TRUE, log.p = FALSE) {
  d     <- .dbase("gompertz", lower.tail = lower.tail, log = log.p, 
                  p = p, shape = shape, rate = rate)
  ret   <- d$ret
  ind   <- d$ind
  p     <- d$p
  shape <- d$shape
  rate  <- d$rate
  
  ret[ind][shape == 0] <- qexp(p[shape == 0], rate = rate[shape == 0])
  sn0 <- shape != 0
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


#' @rdname dgompertz
#' @export
rgompertz <- function(n, shape = 1, rate = 1) {
  r     <- .rbase("gompertz", n = n, shape = shape, rate = rate)
  ret   <- r$ret
  ind   <- r$ind
  shape <- r$shape
  rate  <- r$rate
  
  ret[ind] <- qgompertz(p = runif(sum(ind)), shape = shape, rate = rate)
  ret
}


# == internal helper functions ===============================================


### Standardised procedure for defining density, cumulative
### distribution, hazard and cumulative hazard functions for
### time-to-event distributions

.dbase <- function(dname, lower.tail=TRUE, log=FALSE, ...){
  args <- list(...)
  ## Vectorise all arguments, replicating to length of longest argument
  n <- max(sapply(args, length))
  for (i in seq_along(args)) {
    args[[i]] <- rep(args[[i]], length=n)
  }
  ret <- numeric(n)
  ## Check for parameters out of range, give warning and return NaN
  ## for those
  check.fn <- paste("check.",dname,sep="")
  check.ret <- do.call(check.fn, args[-1])
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

.rbase <- function(dname, n, ...){
  ## Vectorise all arguments, replicating to sample length
  if (length(n) > 1) n <- length(n)
  args <- list(...)
  for (i in seq_along(args)) {
    args[[i]] <- rep(args[[i]], length=n)
  }
  ret <- numeric(n)
  ## Check for parameters out of range, give warning and return NaN
  ## for those
  check.fn <- paste("check.",dname,sep="")
  check.ret <- do.call(check.fn, args)
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
