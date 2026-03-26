
#' Power Calculations for ChiSquared Tests
#' 
#' Compute power of test or determine parameters to obtain target power (same
#' as \code{\link{power.anova.test}}).
#' 
#' Exactly one of the parameters \code{w}, \code{n}, \code{power} or
#' \code{sig.level} must be passed as NULL, and this parameter is determined
#' from the others. Note that the last one has non-NULL default, so \code{NULL}
#' must be explicitly passed, if you want to compute it.
#' 
#' @param n total number of observations.
#' @param w effect size.
#' @param df degree of freedom (depends on the chosen test.
#' @param sig.level Significance level (Type I error probability).
#' @param power Power of test (1 minus Type II error probability).
#' @return Object of class "power.htest", a list of the arguments (including
#' the computed one) augmented with 'method' and 'note' elements.
#' @note \code{\link{uniroot}} is used to solve power equation for unknowns, so
#' you may see errors from it, notably about inability to bracket the root when
#' invalid arguments are given.
#' 
#' @note
#' Based on code by Stephane Champely, and Peter Dalgaard. 
#' 
#' @seealso \code{\link{power.t.test}}
#' @references Cohen, J. (1988) \emph{Statistical power analysis for the
#' behavioral sciences (2nd ed.)} Hillsdale, NJ: Lawrence Erlbaum.
#' @keywords htest
#' @examples
#' 
#' ## Exercise 7.1 P. 249 from Cohen (1988) 
#' powerChisqTest(w=0.289, df=(4-1), n=100, sig.level=0.05)
#' 
#' ## Exercise 7.3 p. 251
#' powerChisqTest(w=0.346, df=(2-1)*(3-1), n=140, sig.level=0.01)
#' 
#' ## Exercise 7.8 p. 270
#' powerChisqTest(w=0.1, df=(5-1)*(6-1), power=0.80, sig.level=0.05)
#' 
 

#' @export
powerChisqTest <- function (n = NULL, w = NULL, df = NULL, sig.level = 0.05, power = NULL) {
  
  if (sum(sapply(list(w, n, df, power, sig.level), is.null)) != 1)
    stop("exactly one of w, n, df, power or sig.level must be NULL")
  if (!is.null(w) && w < 0)
    stop("w must be positive")
  if (!is.null(n) && n < 1)
    stop("number of observations must be at least 1")
  if (!is.null(sig.level) && !is.numeric(sig.level) || any(0 > sig.level | sig.level > 1))
    stop(sQuote("sig.level"), " must be numeric in [0, 1]")
  if (!is.null(power) && !is.numeric(power) || any(0 > power | power > 1))
    stop(sQuote("power"), " must be numeric in [0, 1]")
  p.body <- quote({
    k <- qchisq(sig.level, df = df, lower = FALSE)
    pchisq(k, df = df, ncp = n * w^2, lower = FALSE)
  })
  if (is.null(power))
    power <- eval(p.body)
  else if (is.null(w))
    w <- uniroot(function(w) eval(p.body) - power, c(1e-10, 1e+05))$root
  else if (is.null(n))
    n <- uniroot(function(n) eval(p.body) - power, c(1 + 1e-10, 1e+05))$root
  else if (is.null(sig.level))
    sig.level <- uniroot(function(sig.level) eval(p.body) -
                           power, c(1e-10, 1 - 1e-10))$root
  else stop("internal error")
  
  METHOD <- "Chi squared power calculation"
  NOTE <- "n is the number of observations"
  structure(list(w = w, n = n, df = df, sig.level = sig.level,
                 power = power, method = METHOD, note = NOTE), class = "power.htest")
}



