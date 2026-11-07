

#' Confidence Intervals for a Poisson Rate
#'
#' `poissonCI()` computes confidence intervals for a Poisson event rate
#' using exact and approximate methods.
#'
#' All arguments are vectorized and recycled according to standard R rules.
#'
#' The function assumes \eqn{X \sim \mathrm{Poisson}(n\lambda)}, where
#' \eqn{x} is the observed event count, \eqn{n} the exposure and
#' \eqn{\lambda} the event rate. The point estimate is
#' \eqn{\hat{\lambda} = x/n}.
#'
#' **Exact** (default):
#' The so-called exact interval, based on quantiles of the corresponding
#' chi-square distribution and equivalent to the Garwood (1936) interval.
#' The same interval is returned as `conf.int` by [poisson.test()].
#'
#' **Score**:
#' Obtained by inverting the acceptance region of the Poisson score test.
#'
#' **Wald**:
#' Obtained by inverting the acceptance region of the large-sample normal
#' (Wald) test; the interval is symmetric about \eqn{\hat{\lambda}}.
#'
#' **Byar**:
#' Byar's cube-root normal approximation to the exact interval (Rothman and
#' Boice, 1979). The lower limit is based on \eqn{x}, the upper limit on
#' \eqn{x + 1}; for \eqn{x = 0} the lower limit is 0.
#'
#' Some methods may produce a lower limit below 0. In such cases, the bound
#' is truncated to remain within the valid parameter space
#' \eqn{[0, \infty)}.
#'
#' @section **Which interval should be used?**:
#' The Wald interval has poor coverage for small counts and degenerates to
#' a single point at \eqn{x = 0}. The exact interval guarantees the nominal
#' coverage but is conservative and often unnecessarily wide. The score
#' interval is a good approximate alternative when conservativeness is not
#' required.
#'
#' @param x number of events, a non-negative integer.
#' @param n exposure associated with `x`, such as observation time,
#' person-time or population at risk; a positive number, defaults to 1.
#' @param conf.level confidence level, defaults to 0.95.
#' @param sides a character string specifying the side of the confidence
#' interval, must be one of `"two.sided"` (default), `"left"` or
#' `"right"`. You can specify just the initial letter. `sides` names the
#' side carrying the finite bound: `"left"` reports the lower limit and
#' opens the upper one to `Inf`, `"right"` reports the upper limit and opens
#' the lower one to 0. A one-sided bound at level `conf.level` is the
#' corresponding end of the two-sided interval at level
#' `2 * conf.level - 1`. `"left"` corresponds to the alternative
#' `"greater"` in a hypothesis test.
#' @param method character string specifying which method to use; this can be
#' one out of: `"exact"` (default), `"score"`, `"wald"` and `"byar"`.
#' Several methods may be given at once; all the methods can be asked by
#' `".all"`. Abbreviation of method is accepted. See details.
#'
#' @return If recycling yields a single case, a named numeric vector with
#' elements:
#' \describe{
#'   \item{`est`}{point estimate of the event rate, `x/n`.}
#'   \item{`lci`}{lower confidence interval bound.}
#'   \item{`uci`}{upper confidence interval bound.}
#' }
#'
#' If recycling yields multiple cases, a data frame with one row per case is
#' returned. Its first three columns are `est`, `lci`, and `uci`;
#' the remaining columns contain the recycled argument values.
#'
#' @references
#' Garwood, F. (1936) Fiducial limits for the Poisson distribution.
#' *Biometrika*, **28**, pp. 437-442.
#'
#' Rothman, K. J. and Boice, J. D. Jr. (1979) *Epidemiologic Analysis
#' with a Programmable Calculator*. NIH Publication No. 79-1649. Washington,
#' DC: US Government Printing Office.
#'
#' @seealso [poisson.test()], [poissonDiffCI()], [poissonRatioCI()],
#' [binomCI()]
#'
#' @examples
#'
#' # deaths from horse kicks in 280 Prussian army corps-years
#' count <- 0:4
#' corpsYears <- c(144, 91, 32, 11, 2)
#'
#' x <- sum(count * corpsYears)
#' n <- sum(corpsYears)
#'
#' poissonCI(x, n)
#'
#' # the exact interval corresponds to the one computed by poisson.test
#' poisson.test(x, n)$conf.int
#'
#' # several methods at once
#' poissonCI(x, n, method = c("exact", "score", "wald", "byar"))
#'
#' # a 95% lower confidence bound for the event rate
#' poissonCI(x, n, sides = "left")
#'
#' # SMR for Welsh nickel workers
#' poissonCI(x = 137, n = 24.19893)
#'
#' @concept confidence-interval
#' @concept rate
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
  
  # validated per recycled case, as in .binomCI_engine(); the approximate
  # formulas would otherwise return regular-looking limits for x = 2.5
  checkCount(x, min = 0)

  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) || n <= 0)
    stop(gettextf("'n' must be a positive finite number, got %s",
                  format(n)), call. = FALSE, domain = NA)

  .checkSidedLevel(conf.level, sides)

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

  # Rothman & Boice (1979): the lower limit uses x, the upper limit x + 1,
  # mirroring the chi-square degrees of freedom of the exact interval;
  # at x = 0 the lower limit is 0 (the formula is undefined there)
  byar <- function(k, sgn) k * (1 - 1/(9 * k) + sgn * z/(3 * sqrt(k)))^3

  lci <- if (x == 0) 0 else byar(x, -1) / n
  uci <- byar(x + 1, 1) / n

  return(c(lci = lci, uci = uci))
  
}
