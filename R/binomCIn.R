
#' Sample Size for a Desired Width of a Binomial Confidence Interval
#'
#' `binomCIn()` computes the required sample size to obtain a binomial
#' confidence interval of a specified width, as calculated by `binomCI()`.
#' The function uses [uniroot()] to numerically solve for the
#' corresponding sample size.
#'
#' **Required Samplesize** (by `binomCIn()`): \cr
#' The required sample size for a given confidence interval width depends
#' on the assumed population proportion. Since this proportion is often
#' unknown at the planning stage of a study, a conservative approach is to
#' use the worst-case scenario of \eqn{p = 0.5}, which maximizes the variance
#' and therefore yields the largest required sample size. If a more 
#' accurate estimate of the population proportion is available,
#' it can be used to obtain a smaller required sample size for the same
#' level of precision.
#'
#' The root search evaluates the interval at \eqn{x = p \cdot n} for
#' continuous \eqn{n}, so only those methods can be inverted whose limits
#' are smooth functions of the count. The methods `"mid-p"`, `"blaker"`,
#' `"witting"` and `"likelihood"` are defined through the discrete binomial
#' distribution (and `"witting"` is randomized on top of that); they are
#' rejected with an error instead of returning a silently meaningless root.
#'
#' The returned sample size is not rounded, round it up to get a feasible
#' number of observations.
#'  
#' @param p probability for success, defaults to `0.5` as worst case.
#' @param width the width of the confidence interval.
#' @param interval a vector containing the end-points of the interval to be
#' searched for the root. The defaults are set to `c(1, 100000)`.
#'
#' @return `binomCIn` returns a single numeric value giving the required sample size.
#' 
#' @examples
#' 
#' binomCIn(p=0.1, width=0.05, method="pratt")
#' 
#' # round up to get the number of observations to plan for
#' ceiling(binomCIn(width=0.1))
#' 

#' @rdname binomCI

#' @family ci.proportion  
#' @concept confidence-interval  
#' @concept proportion  
#' @concept sample-size
#'
#'
#' @export
binomCIn <- function(p=0.5, width, interval=c(1, 1e5), 
                     conf.level=0.95, sides=c("two.sided", "left", "right"),
                     method="wilson") {
  
  sides <- match.arg(sides)
  method <- match.arg(method, eval(formals(binomCI)$method))
  
  # the root search evaluates the interval at non-integer x = p * n, which
  # only the methods with closed form limits tolerate
  if (method %in% c("mid-p", "blaker", "witting", "likelihood"))
    stop(gettextf(paste("method '%s' is a function of integer counts and",
                        "cannot be inverted for a sample size"),
                  method), call. = FALSE)
  
  if (!isTRUE(is.numeric(p) && length(p) == 1L && p > 0 && p < 1))
    stop("'p' must be a single probability in (0, 1)", call. = FALSE)
  
  if (!isTRUE(is.numeric(width) && length(width) == 1L && width > 0 && width < 1))
    stop("'width' must be a single value in (0, 1)", call. = FALSE)
  
  if (!isTRUE(is.numeric(interval) && length(interval) == 2L &&
              all(is.finite(interval)) && interval[1L] > 0 &&
              interval[1L] < interval[2L]))
    stop("'interval' must be two increasing positive values", call. = FALSE)
  
  checkConfLevel(conf.level)
  if (is.na(conf.level))
    stop("'conf.level' must not be NA when solving for a sample size",
         call. = FALSE)
  
  alpha <- .sidesAlpha(conf.level, sides)
  
  # go directly for the bare limits, binomCI() would (rightly) reject the
  # non-integer x the root search needs
  f <- function(n) {
    ci <- applySides(unname(.binomCI_bounds(p * n, n, alpha, method)[c("lci", "uci")]),
                     sides, lo = 0, hi = 1)
    ci[["uci"]] - ci[["lci"]] - width
  }
  
  # at very small n some limits are undefined for non-integer x = p * n
  # (pratt needs n * (1 - p) > 1, for instance). The interval is clamped to
  # [0, 1] there anyway, so the lower end can be moved up without losing the
  # root
  lo <- interval[1L]
  while (!isTRUE(is.finite(suppressWarnings(f(lo)))) && lo < interval[2L])
    lo <- min(lo * 2, interval[2L])
  
  if (!isTRUE(is.finite(suppressWarnings(f(lo)))))
    stop(gettextf("the width of the '%s' interval is not computable on 'interval'",
                  method), call. = FALSE)
  
  uniroot(f, interval = c(lo, interval[2L]))$root
  
}

