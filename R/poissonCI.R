

#' Confidence Interval for a Poisson Rate
#'
#' Estimates a Poisson event rate and calculates an exact or approximate
#' confidence interval.
#'
#' @param x non-negative integer event count or vector of counts
#' @param n positive exposure associated with `x`, such as observation time,
#'   person-time, or population at risk; may be a vector and defaults to 1
#' @param conf.level numeric confidence level between 0 and 1; defaults to 0.95
#' @param sides type of confidence interval: `"two.sided"`, `"left"`, or
#'   `"right"`; may be abbreviated
#' @param method method used to calculate the confidence interval: `"exact"`,
#'   `"score"`, `"wald"`, or `"byar"`; may be abbreviated and may contain
#'   several methods; defaults to `"exact"`
#'
#' @return If all arguments identify a single result, a named numeric vector
#'   with elements:
#'   \describe{
#'     \item{`est`}{estimated event rate}
#'     \item{`lci`}{lower confidence bound}
#'     \item{`uci`}{upper confidence bound}
#'   }
#'   Otherwise, a `data.frame` containing these three columns and the recycled
#'   argument values that identify each result.
#'
#' @details
#' The function assumes
#' \deqn{X \sim \mathrm{Poisson}(n\lambda),}
#' where `x` is the observed event count, `n` is the exposure, and
#' \eqn{\lambda} is the event rate. The point estimate is
#' \eqn{\hat{\lambda} = x/n}.
#'
#' The available confidence-interval methods are:
#' \describe{
#'   \item{`"exact"`}{the exact Poisson interval calculated by
#'     [stats::poisson.test()], equivalent to the Garwood interval}
#'   \item{`"score"`}{the interval obtained by inverting the Poisson score
#'     test}
#'   \item{`"wald"`}{the normal-approximation interval centred at
#'     \eqn{\hat{\lambda}}}
#'   \item{`"byar"`}{Byar's cube-root normal approximation}
#' }
#'
#' The lower bound is restricted to the parameter space \eqn{[0, \infty)}.
#' For `sides = "left"`, the function returns a lower one-sided confidence
#' bound and sets `uci` to `Inf`. This corresponds to the alternative
#' `"greater"` in a hypothesis test. For `sides = "right"`, it returns an
#' upper one-sided confidence bound and sets `lci` to 0.
#'
#' Compatible vector lengths are recycled. Supplying several methods therefore
#' provides the corresponding intervals in a single `data.frame`.
#'
#' @references
#' Garwood, F. (1936). Fiducial limits for the Poisson distribution.
#' \emph{Biometrika}, \bold{28}, 437--442.
#'
#' Rothman, K. J. and Boice, J. D. Jr. (1979). \emph{Epidemiologic Analysis
#' with a Programmable Calculator}. NIH Publication No. 79-1649. Washington,
#' DC: US Government Printing Office.
#'
#' @seealso [stats::poisson.test()]
#'
#' @examples
#' # Deaths from horse kicks in 280 Prussian army corps-years
#' count <- 0:4
#' corpsYears <- c(144, 91, 32, 11, 2)
#'
#' x <- sum(count * corpsYears)
#' n <- sum(corpsYears)
#'
#' poissonCI(x, n)
#' poissonCI(x, n, method = c("exact", "score", "wald", "byar"))
#'
#' # A 95% lower confidence bound for the event rate
#' poissonCI(x, n, sides = "left")
#'
#' # SMR for Welsh nickel workers
#' poissonCI(x = 137, n = 24.19893)
#'
#' @concept confidence-interval
#'
#' @export
poissonCI <- function(x, n = 1, conf.level = 0.95, 
                      sides = c("two.sided","left","right"),
                      method = c("exact","score", "wald","byar")) {

  sides <- match.arg(sides)
  
  if (missing(method)) {
    # if not provided take the first method instead of all (!)
    method <- eval(formals(sys.function())$method)[1]
    
  } else {
    # resolve methods cleanly, allowing an ".all" hidden option for method
    method <- .resolveMethod(method, several.ok = TRUE)
  }
  
  res <- .recycleApply(.poissonCI_engine,
                       x = x,
                       n = n,
                       conf.level = conf.level,
                       sides = sides,
                       method = method
  )
  
  if(length(res) == 1)
    out <- res[[1]]
  else{
    out <- as.data.frame(attr(res, "recycle"))
    out <- data.frame(do.call(rbind, res), out)
  }
  
  return(out)
  
}  


# ==  internal helper functions  ===========================================


#' @keywords internal
.poissonCI_engine <- function(x, n, conf.level, sides, method, stdEst = NULL){

  # ref:  http://www.ijmo.org/papers/189-S083.pdf but wacklig!!!
  # http://www.math.montana.edu/~rjboik/classes/502/ci.pdf
  # http://www.ine.pt/revstat/pdf/rs120203.pdf
  # http://www.pvamu.edu/include/Math/AAM/AAM%20Vol%206,%20Issue%201%20(June%202011)/06_%20Kibria_AAM_R308_BK_090110_Vol_6_Issue_1.pdf
  
  # see also:   pois.conf.int {epitools}
  
  alpha <- 1 - conf.level
  
  # For a one-sided interval, the reported bound is one side of a
  # (narrower, two-sided) interval computed at a DOUBLED alpha: the
  # two-sided interval's tail probability alpha_used/2 must equal the
  # target one-sided alpha, i.e. alpha_used = 2 * alpha. Verified against
  # poisson.test(..., alternative = "greater")'s native one-sided bound.
  if (sides != "two.sided")
    alpha <- alpha * 2
  
  CI <- switch( method
                , "exact" = { .poissonCI.exact(x, n, alpha) }
                , "wald" =  { .poissonCI.wald(x, n, alpha) }
                , "score" = { .poissonCI.score(x, n, alpha) }
                , "byar" =  { .poissonCI.byar(x, n, alpha) }
                , stop("Unknown method.")
      # agresti-coull is the same as score
      # garwood is the same as exact, check that!!
  )
  
  # this is the default lambda estimator
  est <- x/n
  
  # dot not return ci bounds outside [0, Inf]
  ci <- c( est = est, 
           lci = max(0, unname(CI["lci"])), 
           uci = unname(CI["uci"]) )    # no limits on the right side
  
  if(sides=="left")
    ci[3] <- Inf
  else if(sides=="right")
    ci[2] <- 0
  
  return(ci)
  
}



#' @keywords internal
.poissonCI.exact <- function(x, n, alpha) {
    
    return(setNamesX(
        poisson.test(x, n, conf.level = 1-alpha)$conf.int, 
      c("lci","uci")))
}


#' @keywords internal
.poissonCI.score <- function(x, n, alpha) {
  
  z <- qnorm(1 - alpha/2)

  term1 <- (x + z^2/2)/n
  term2 <- z * n^-0.5 * sqrt(x/n + z^2/(4*n))

  return(setNamesX(term1 - c(1,-1) * term2, c("lci","uci")))
  
}


#' @keywords internal
.poissonCI.wald <- function(x, n, alpha) {
  
  z <- qnorm(1-alpha/2)
  lambda <- x/n

  term2 <- z*sqrt(lambda/n)
  return(setNamesX(lambda - c(1,-1) * term2, c("lci","uci")))
  
}

#' @keywords internal
.poissonCI.byar <- function(x, n, alpha) {
  
  z <- qnorm(1-alpha/2)

  xcc <- x + 0.5
  zz  <- (z/3) * sqrt(1/xcc)
  
  lci <- (xcc * (1 - 1/(9 * xcc) - zz)^3)/n
  uci <- (xcc * (1 - 1/(9 * xcc) + zz)^3)/n
  
  return(c(lci=lci, uci=uci))
  
}
