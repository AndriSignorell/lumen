
#' Confidence Intervals for the Number of Successes in a Finite Population
#'
#' `hyperCI()` computes confidence intervals for the number \eqn{M} of
#' successes in a finite population of size \eqn{N}, when \eqn{x} successes
#' are observed in a sample of size \eqn{n} drawn without replacement, i.e.
#' \eqn{X \sim \mathrm{Hyper}(M, N - M, n)}. It is the finite population
#' counterpart of [binomCI()]. The limits are counts; divide them by `N` for
#' the population proportion.
#'
#' All arguments are vectorized and recycled according to standard R rules.
#'
#' The sample already fixes \eqn{x \le M \le N - n + x}. All limits are kept
#' within this range, and a one-sided interval opens its free side to it.
#'
#' **Wald**:
#' The Wald interval for the proportion with the finite population
#' correction, \eqn{\hat p \pm z \sqrt{f \hat p (1 - \hat p) / n}} with
#' \eqn{\hat p = x/n} and \eqn{f = (N - n)/(N - 1)} (Cochran 1977),
#' multiplied by \eqn{N}.
#'
#' **Wilson** (default):
#' The Wilson score interval of [binomCI()] with the effective sample size
#' \eqn{n / f}, which carries the finite population correction into the score
#' variance, multiplied by \eqn{N}.
#'
#' The two asymptotic intervals are rounded outwards to integers. For a census
#' (\eqn{n = N}) they return \eqn{M = x}.
#'
#' **Clopper-Pearson**:
#' The exact interval obtained by inverting two one-sided hypergeometric tests
#' (Konijn 1973): the lower limit is the smallest \eqn{M} with
#' \eqn{P(X \ge x \mid M) > \alpha/2}, the upper limit the largest \eqn{M}
#' with \eqn{P(X \le x \mid M) > \alpha/2}. It guarantees the coverage but is
#' conservative.
#'
#' **Mid-p**:
#' As Clopper-Pearson, with half the probability of the observed count in
#' each tail. Not guaranteed to reach the level, but closer to it on average.
#'
#' **Blaker**:
#' The exact interval of Blaker (2000), obtained by inverting the test that
#' uses the smaller tail probability as statistic. It is contained in the
#' Clopper-Pearson interval and keeps the level; its acceptance region can in
#' rare cases have gaps, the interval spans the smallest and largest accepted
#' value. One-sided, the Clopper-Pearson bound is returned: an end of the
#' two-sided Blaker interval at level `2 * conf.level - 1` misses the level.
#'
#' **Wang**:
#' The admissible exact interval of Wang (2015). Starting from the
#' Clopper-Pearson interval, the limits are shrunk pairwise
#' (\eqn{U_x = N - L_{n-x}}) from the middle of the sample space outwards, each
#' as far as the coverage permits. The resulting family is monotone and
#' symmetric, and no limit can be moved inwards, the other intervals held
#' fixed, without the coverage falling below `conf.level`. It is never wider
#' than Clopper-Pearson. One-sided, the Clopper-Pearson bound is already the
#' smallest exact bound and is returned instead of an end of the two-sided
#' Wang interval at level `2 * conf.level - 1`, which would miss the level.
#' The computation proceeds from \eqn{n/2} towards \eqn{x}; it takes below a
#' second for \eqn{n = 5000}, \eqn{N = 10^6}.
#'
#' @param x number of successes in the sample, an integer between 0 and `n`.
#' @param n sample size, a positive integer not larger than `N`.
#' @param N population size, a positive integer.
#' @param conf.level confidence level, defaults to 0.95. With `NA` only the
#' point estimate is returned.
#' @param sides a character string specifying the side of the confidence
#' interval, must be one of `"two.sided"` (default), `"left"` or
#' `"right"`. You can specify just the initial letter. `sides` names the
#' side carrying the finite bound: `"left"` reports the lower limit and
#' opens the upper one to `N - n + x`, `"right"` reports the upper limit and
#' opens the lower one to `x`. A one-sided bound at level `conf.level` is the
#' corresponding end of the two-sided interval at level
#' `2 * conf.level - 1`, and therefore requires `conf.level > 0.5`. The
#' exceptions are `"blaker"` and `"wang"`, which are calibrated on the
#' two-sided coverage only; their one-sided bound is the Clopper-Pearson
#' bound (see details).
#' @param method character string specifying which method to use; this can be
#' one out of: `"wilson"` (default), `"wald"`, `"clopper-pearson"`,
#' `"mid-p"`, `"blaker"` and `"wang"`. All the methods can be asked by
#' `".all"`. Abbreviation of method is accepted. See details.
#'
#' @return If recycling yields a single case, a named numeric vector with
#' elements:
#' \describe{
#'   \item{`est`}{point estimate of the number of successes in the
#'     population, the unbiased \eqn{N x / n}.}
#'   \item{`lci`}{lower confidence interval bound, an integer.}
#'   \item{`uci`}{upper confidence interval bound, an integer.}
#' }
#'
#' If recycling yields multiple cases, a data frame with one row per case is
#' returned. Its first three columns are `est`, `lci`, and `uci`;
#' the remaining columns contain the recycled argument values.
#'
#' With `conf.level = NA` no interval is computed: the point estimate is
#' returned as an unnamed scalar, or as the single column `est` of the
#' data frame.
#'
#' @references
#' Blaker, H. (2000) Confidence curves and improved exact confidence intervals
#' for discrete distributions, *Canadian Journal of Statistics* 28 (4),
#' 783-798
#'
#' Cochran, W. G. (1977) *Sampling Techniques*, 3rd ed. New York: Wiley.
#'
#' Konijn, H. S. (1973) *Statistical Theory of Sample Survey Design and
#' Analysis*. Amsterdam: North-Holland.
#'
#' Wang, W. (2015) Exact optimal confidence intervals for hypergeometric
#' parameters, *Journal of the American Statistical Association* 110 (512),
#' 1491-1499, \doi{10.1080/01621459.2014.966191}
#'
#' Wilson, E. B. (1927) Probable inference, the law of succession, and
#' statistical inference, *Journal of the American Statistical
#' Association* 22, 209-212.
#'
#' @seealso [binomCI()] for sampling with replacement or an infinite
#' population, [phyper()]
#'
#' @examples
#' # audit: 10 faulty invoices in a sample of 50 out of 2000
#' hyperCI(x = 10, n = 50, N = 2000, method = ".all")
#'
#' # the same as proportions
#' hyperCI(x = 10, n = 50, N = 2000, method = "wang")[c("lci", "uci")] / 2000
#'
#' # the finite population correction matters once n is a sizeable
#' # fraction of N, and vanishes for large N
#' hyperCI(x = 10, n = 50, N = c(60, 200, 1e6), method = "clopper-pearson")
#' binomCI(x = 10, n = 50, method = "clopper-pearson")
#'
#' # an upper bound for the number of faulty items after a clean sample
#' hyperCI(x = 0, n = 50, N = 2000, sides = "right", method = "clopper-pearson")
#'
#' @family ci.proportion
#' @concept confidence-interval
#' @concept proportion
#'
#'
#' @export
hyperCI <- function(x, n, N,
                    conf.level = 0.95,
                    sides = c("two.sided", "left", "right"),
                    method = c("wilson", "wald", "clopper-pearson",
                               "mid-p", "blaker", "wang")) {

  # conf.level is validated in the engine, where it arrives as a single
  # recycled value; 'sides' is resolved here because applySides() expects
  # a matched value
  sides <- match.arg(sides)

  if (missing(method)) {
    # if not provided take the first method instead of all (!)
    method <- eval(formals(sys.function())$method)[1]

  } else {
    # resolve methods cleanly, allowing an ".all" hidden option for method
    method <- .resolveMethod(method, several.ok = TRUE)
  }

  res <- .recycleApply(.hyperCI_engine,
                       x = x,
                       n = n,
                       N = N,
                       conf.level = conf.level,
                       sides = sides,
                       method = method
                       )

  if(length(res) == 1)
    out <- res[[1]]
  else{
    ci <- do.call(rbind, res)
    # with conf.level = NA the engine returns the bare point estimate
    if(is.null(colnames(ci)))
      colnames(ci) <- "est"
    out <- data.frame(ci, as.data.frame(attr(res, "recycle")))
  }

  return(out)

}


# ==  internal helper functions  ===========================================

#' @keywords internal
.hyperCI_engine <- function(x, n, N, conf.level, sides, method){

  checkCount(N, min = 1)
  checkCount(n, min = 1)
  checkCount(x, min = 0)

  if (n > N)
    stop(gettextf("'n' must not be larger than 'N', got n = %g and N = %g",
                  n, N), call. = FALSE, domain = NA)

  if (x > n)
    stop(gettextf("'x' must not be larger than 'n', got x = %g and n = %g",
                  x, n), call. = FALSE, domain = NA)

  checkConfLevel(conf.level)

  if (sides != "two.sided" && !is.na(conf.level) && conf.level <= 0.5)
    stop(gettextf(
      "a one-sided interval needs 'conf.level' above 0.5, not %g",
      conf.level), domain = NA)

  est <- N * x / n

  if (is.na(conf.level))
    return(unname(est))

  alpha <- .sidesAlpha(conf.level, sides)

  CI <- .hyperCI_bounds(x, n, N, alpha, method, sides)

  # the sample fixes x <= M <= N - n + x; clamping to this range and opening
  # the free side happen in one place (design_rules 8.1.4)
  c(est = unname(est),
    applySides(c(CI[["lci"]], CI[["uci"]]), sides, lo = x, hi = N - n + x))

}


#' @keywords internal
# the naked interval for M: no validation, no clamping, no point estimate.
# 'sides' is only needed for the methods whose one-sided bound is not an end
# of the two-sided interval at the doubled alpha (blaker, wang).
.hyperCI_bounds <- function(x, n, N, alpha, method, sides = "two.sided") {

  # blaker and wang are calibrated on the two-sided coverage only: one end of
  # their two-sided interval at the doubled alpha is not a one-sided bound at
  # conf.level and undercovers. One-sided, the Clopper-Pearson bound is the
  # smallest exact bound and is used instead.
  if (sides != "two.sided" && method %in% c("blaker", "wang"))
    method <- "clopper-pearson"

  switch( method
          , "wald" =              { .hyperCI.asymp(x, n, N, alpha, "wald") }
          , "wilson" =            { .hyperCI.asymp(x, n, N, alpha, "wilson") }
          , "clopper-pearson" =   { .hyperCI.exact(x, n, N, alpha) }
          , "mid-p" =             { .hyperCI.exact(x, n, N, alpha, midp = TRUE) }
          , "blaker" =            { .hyperCI.blaker(x, n, N, alpha) }
          , "wang" =              { .hyperCI.wang(x, n, N, alpha) }
          , stop(gettextf("Unknown method '%s'.", method))
  )

}


#' @keywords internal
# wald and wilson from binomCI with the effective sample size n / f, where
# f = (N - n) / (N - 1) is the finite population correction (x/n * n/f = x/f
# successes), scaled to counts and rounded outwards
.hyperCI.asymp <- function(x, n, N, alpha, method) {

  f <- if (N > 1) (N - n) / (N - 1) else 0

  # census: the sample is the population
  if (f == 0)
    return(c(lci = x, uci = x))

  ci <- .binomCI_bounds(x / f, n / f, alpha, method)

  # the tolerance keeps an integer limit from being rounded outwards by the
  # last bits of the floating point product
  c(lci = floor(N * ci[["lci"]] + 1e-8),
    uci = ceiling(N * ci[["uci"]] - 1e-8))

}


#' @keywords internal
# smallest integer in [lo, hi] for which the nondecreasing predicate 'ok'
# holds; 'ok(hi)' must be TRUE
.hyperFirst <- function(ok, lo, hi) {
  while (lo < hi) {
    mid <- lo + (hi - lo) %/% 2
    if (ok(mid)) hi <- mid else lo <- mid + 1
  }
  lo
}


#' @keywords internal
# Clopper-Pearson type (Konijn 1973) and mid-p limits. The lower limit is the
# smallest M whose upper tail at x exceeds alpha/2; the tail is nondecreasing
# in M and reaches 1 (mid-p: at least 1/2) at M = N - n + x. The upper limit
# follows from the symmetry n - X ~ Hyper(N - M, M, n).
.hyperCI.exact <- function(x, n, N, alpha, midp = FALSE) {

  lower <- function(x)
    .hyperFirst(function(M)
      phyper(x - 1, M, N - M, n, lower.tail = FALSE) -
        midp * 0.5 * dhyper(x, M, N - M, n) > alpha / 2,
      x, N - n + x)

  c(lci = lower(x), uci = N - lower(n - x))

}


#' @keywords internal
# Blaker's acceptability of x under M: the probability of all counts whose
# smaller tail is not larger than that of x (relative tolerance for ties as
# in exactci)
.hyperBlakerAccept <- function(x, n, N, M) {

  y  <- max(0, n - N + M):min(n, M)
  d  <- dhyper(y, M, N - M, n)
  gm <- pmin(cumsum(d), rev(cumsum(rev(d))))   # min(P(X <= y), P(X >= y))

  sum(d[gm <= gm[y == x] * (1 + 1e-7)])

}


#' @keywords internal
# the Blaker interval is contained in the Clopper-Pearson interval, so the
# lower limit is found by scanning upwards from the CP limit to the first
# accepted value; the upper limit by symmetry
.hyperCI.blaker <- function(x, n, N, alpha) {

  cp <- .hyperCI.exact(x, n, N, alpha)

  lower <- function(x, M) {
    while (M < N - n + x && .hyperBlakerAccept(x, n, N, M) <= alpha)
      M <- M + 1
    M
  }

  c(lci = lower(x, cp[["lci"]]),
    uci = N - lower(n - x, N - cp[["uci"]]))

}


#' @keywords internal
# Wang (2015): admissible exact interval, see src/wangHyperCI.cpp. Two-sided
# only; .hyperCI_bounds() hands one-sided requests to clopper-pearson.
.hyperCI.wang <- function(x, n, N, alpha) {
  setNames(as.numeric(.wangHyperCI(x, n, N, alpha)), c("lci", "uci"))
}
