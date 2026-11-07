
# Exact (Garwood) limits for a single Poisson rate at one-sided level a.
# Identical to poissonCI(method = "exact"), inlined here so that the recycling
# loop does not pay for a full htest-style call per row.
.poissonExactLo <- function(x, n, a) if(x == 0) 0 else stats::qchisq(a, 2 * x) / (2 * n)
.poissonExactHi <- function(x, n, a) stats::qchisq(1 - a, 2 * (x + 1)) / (2 * n)


# Root of a mid-p equation h(p) = 0, searched on the log-odds scale and
# returned there: the rate ratio is f * exp(eta), so an absolute tolerance on
# eta is a relative tolerance on the limit, also for very unbalanced counts
# where p is close to 0 or 1. At extreme one-sided levels the equation can
# have no root at all, in which case the boundary nearest to solving it is
# the honest answer and maps to a limit of 0 (-Inf) or Inf (Inf).
.midpSolve <- function(h) {
  g <- function(eta) h(stats::plogis(eta))
  lo <- stats::qlogis(.Machine$double.eps)
  hi <- -lo
  flo <- g(lo)
  fhi <- g(hi)
  if(!is.finite(flo) || !is.finite(fhi) || flo * fhi > 0)
    return(if(abs(flo) <= abs(fhi)) -Inf else Inf)
  stats::uniroot(g, c(lo, hi), f.lower = flo, f.upper = fhi,
                 tol = 1e-10)$root
}


# Shared recycling wrapper: returns a named vector for a single case and a
# data frame with one row per case otherwise, as in binomDiffCI/binomRatioCI.
.poissonTwoSampleCI <- function(x1, n1, x2, n2, conf.level, sides, method, fun) {

  args <- list(x1 = x1, n1 = n1, x2 = x2, n2 = n2, conf.level = conf.level)

  if(!all(vapply(args, is.numeric, logical(1))))
    stop("'x1', 'n1', 'x2', 'n2' and 'conf.level' must be numeric")
  if(any(lengths(args) == 0L))
    stop("arguments must not have length zero")
  if(!all(vapply(args, function(z) all(is.finite(z)), logical(1))))
    stop("'x1', 'n1', 'x2', 'n2' and 'conf.level' must be finite")
  if(any(c(x1, x2) < 0) || any(c(x1, x2) != round(c(x1, x2))))
    stop("'x1' and 'x2' must be non-negative integers")
  if(any(c(n1, n2) <= 0))
    stop("'n1' and 'n2' must be positive")
  if(any(conf.level <= 0) || any(conf.level >= 1))
    stop("'conf.level' must lie in (0, 1)")

  # recycle only to whole multiples, so incompatible lengths are an error
  # rather than a silent truncation
  len <- max(lengths(args))
  if(any(len %% lengths(args) != 0L))
    stop("arguments cannot be recycled to a common length")
  args <- lapply(args, rep, length.out = len)

  res <- t(vapply(seq_len(len),
                  function(i) fun(args$x1[i], args$n1[i], args$x2[i], args$n2[i],
                                  args$conf.level[i], sides, method),
                  numeric(3)))
  colnames(res) <- c("est", "lci", "uci")

  if(len == 1L)
    return(res[1L, ])

  data.frame(res, as.data.frame(args))
}


#' Confidence Intervals for the Difference Between Two Poisson Rates
#'
#' `poissonDiffCI()` computes confidence intervals for the difference between
#' two independent Poisson event rates.
#'
#' The numeric arguments are recycled to a common length only when their
#' lengths are compatible whole multiples; incompatible lengths produce an
#' error. `sides` and `method` must each identify a single choice.
#'
#' The function assumes two independent counts
#' \eqn{X_i \sim \mathrm{Poisson}(n_i\lambda_i),\ i = 1, 2}. The parameter of
#' interest is the rate difference \eqn{\Delta = \lambda_1 - \lambda_2},
#' estimated by \eqn{\hat{\Delta} = x_1/n_1 - x_2/n_2}.
#'
#' **MOVER** (default):
#' The method of variance estimates recovery (Zou and Donner, 2008), which
#' combines the separate exact Garwood limits of the two rates (see
#' [poissonCI()]).
#'
#' **Wald**:
#' Obtained by inverting the acceptance region of the large-sample normal
#' (Wald) test, based on the standard error
#' \eqn{\sqrt{x_1/n_1^2 + x_2/n_2^2}}. The interval has zero width when both
#' counts are zero.
#'
#' @section **Which interval should be used?**:
#' The Wald interval inherits the poor small-count behaviour of the
#' one-sample Wald interval. The MOVER interval builds on the exact limits of
#' the two rates and is the better choice for small or unbalanced counts
#' (Zou and Donner, 2008).
#'
#' @param x1 number of events in the first sample, a non-negative integer.
#' @param n1 exposure associated with `x1`, such as observation time,
#' person-time or population at risk; a positive number, defaults to 1.
#' @param x2 number of events in the second sample, a non-negative integer.
#' @param n2 exposure associated with `x2`; a positive number, defaults to 1.
#' @param conf.level confidence level, defaults to 0.95.
#' @param sides a character string specifying the side of the confidence
#' interval, must be one of `"two.sided"` (default), `"left"` or
#' `"right"`. You can specify just the initial letter. `sides` names the
#' side carrying the finite bound: `"left"` reports the lower limit and
#' opens the upper one to `Inf`, `"right"` reports the upper limit and opens
#' the lower one to `-Inf`. A one-sided bound puts the whole
#' `1 - conf.level` into its single tail. One-sided `"mover"` intervals
#' require `conf.level >= 0.5`: below, the Garwood limits cross the point
#' estimates and the combined bound would no longer be monotone in the level.
#' @param method character string specifying which method to use; this can be
#' one out of: `"mover"` (default) and `"wald"`. Abbreviation of method is
#' accepted. See details.
#'
#' @return If recycling yields a single case, a named numeric vector with
#' elements:
#' \describe{
#'   \item{`est`}{point estimate of the rate difference.}
#'   \item{`lci`}{lower confidence interval bound.}
#'   \item{`uci`}{upper confidence interval bound.}
#' }
#'
#' If recycling yields multiple cases, a data frame with one row per case is
#' returned. Its first three columns are `est`, `lci`, and `uci`;
#' the remaining columns contain the recycled values of `x1`, `n1`, `x2`,
#' `n2`, and `conf.level`.
#'
#' @references
#' Zou, G. Y. and Donner, A. (2008) Construction of confidence limits about
#' effect measures: a general approach. *Statistics in Medicine*,
#' **27**(10), pp. 1693-1702.
#'
#' @seealso [poissonCI()], [poissonRatioCI()], [binomDiffCI()]
#'
#' @examples
#'
#' # 15 events in 100 person-years compared with
#' # 6 events in 120 person-years
#' poissonDiffCI(15, 100, 6, 120)
#' poissonDiffCI(15, 100, 6, 120, method = "wald")
#'
#' # a 95% lower confidence bound for the rate difference
#' poissonDiffCI(15, 100, 6, 120, sides = "left")
#'
#' # recycling returns one row per comparison
#' poissonDiffCI(x1 = c(15, 20), n1 = 100, x2 = 6, n2 = 120)
#'
#' @concept confidence-interval
#' @concept rate
#'
#' @export
poissonDiffCI <- function(x1, n1 = 1, x2, n2 = 1, conf.level = 0.95,
                          sides = c("two.sided", "left", "right"),
                          method = c("mover", "wald")) {

  sides <- match.arg(sides)
  method <- match.arg(method)

  .poissonTwoSampleCI(x1, n1, x2, n2, conf.level, sides, method,
                      fun = function(x1, n1, x2, n2, conf.level, sides, method) {

    # below 0.5 the one-sided Garwood limits cross the point estimates, and
    # squaring the distances drops that sign: the bound would move outwards
    # as the level falls
    if(method == "mover" && sides != "two.sided" && conf.level < 0.5)
      stop(gettextf(
        "one-sided 'mover' intervals need 'conf.level' of at least 0.5, not %g",
        conf.level), call. = FALSE, domain = NA)

    a <- if(sides == "two.sided") (1 - conf.level) / 2 else 1 - conf.level
    r1 <- x1 / n1
    r2 <- x2 / n2
    est <- r1 - r2

    if(method == "wald") {
      h <- stats::qnorm(1 - a) * sqrt(x1 / n1^2 + x2 / n2^2)
      lci <- est - h
      uci <- est + h
    } else {
      l1 <- .poissonExactLo(x1, n1, a); u1 <- .poissonExactHi(x1, n1, a)
      l2 <- .poissonExactLo(x2, n2, a); u2 <- .poissonExactHi(x2, n2, a)
      lci <- est - sqrt((r1 - l1)^2 + (u2 - r2)^2)
      uci <- est + sqrt((u1 - r1)^2 + (r2 - l2)^2)
    }

    if(sides == "left")  uci <- Inf
    if(sides == "right") lci <- -Inf

    c(est, lci, uci)
  })
}


#' Confidence Intervals for the Ratio of Two Poisson Rates
#'
#' `poissonRatioCI()` computes confidence intervals for the ratio of two
#' independent Poisson event rates.
#'
#' The numeric arguments are recycled to a common length only when their
#' lengths are compatible whole multiples; incompatible lengths produce an
#' error. `sides` and `method` must each identify a single choice.
#'
#' The function assumes two independent counts
#' \eqn{X_i \sim \mathrm{Poisson}(n_i\lambda_i),\ i = 1, 2}. The parameter of
#' interest is the rate ratio \eqn{\theta = \lambda_1/\lambda_2}, estimated
#' by \eqn{\hat{\theta} = (x_1/n_1)/(x_2/n_2)}.
#'
#' **Exact** (default):
#' The exact conditional interval, obtained by conditioning on
#' \eqn{x_1 + x_2} and transforming a Clopper-Pearson interval for the
#' resulting binomial probability. This is the construction used by
#' [poisson.test()] for two samples.
#'
#' **Mid-p**:
#' The corresponding conditional mid-p interval, designed to reduce the
#' conservatism of the exact interval. It is generally shorter but does not
#' guarantee the nominal coverage.
#'
#' **Wald-log**:
#' The asymptotic Wald interval, symmetric on the logarithmic scale. It
#' cannot be calculated when either count is zero and is then returned as
#' \eqn{[0, \infty]}.
#'
#' If both counts are zero, the rate ratio is undefined; the function returns
#' `NA` for `est` and \eqn{[0, \infty]} for the confidence interval for every
#' method. If only \eqn{x_2} is zero, the point estimate is `Inf`.
#'
#' @section **Which interval should be used?**:
#' The exact interval guarantees the nominal coverage but is conservative.
#' The mid-p interval is a common compromise with coverage close to the
#' nominal level (Graham et al., 2003). The log-Wald interval requires
#' reasonably large counts in both samples and is not available when either
#' count is zero.
#'
#' @inheritParams poissonDiffCI
#' @param sides a character string specifying the side of the confidence
#' interval, must be one of `"two.sided"` (default), `"left"` or
#' `"right"`. You can specify just the initial letter. `sides` names the
#' side carrying the finite bound: `"left"` reports the lower limit and
#' opens the upper one to `Inf`, `"right"` reports the upper limit and opens
#' the lower one to 0, the lower limit of the parameter space. A one-sided
#' bound puts the whole `1 - conf.level` into its single tail.
#' @param method character string specifying which method to use; this can be
#' one out of: `"exact"` (default), `"mid-p"` and `"wald-log"`.
#' Abbreviation of method is accepted. See details.
#'
#' @return If recycling yields a single case, a named numeric vector with
#' elements:
#' \describe{
#'   \item{`est`}{point estimate of the rate ratio.}
#'   \item{`lci`}{lower confidence interval bound.}
#'   \item{`uci`}{upper confidence interval bound.}
#' }
#'
#' If recycling yields multiple cases, a data frame with one row per case is
#' returned. Its first three columns are `est`, `lci`, and `uci`;
#' the remaining columns contain the recycled values of `x1`, `n1`, `x2`,
#' `n2`, and `conf.level`.
#'
#' @references
#' Graham, P. L., Mengersen, K. and Morton, A. P. (2003) Confidence limits
#' for the ratio of two rates based on likelihood scores. *Statistics in
#' Medicine*, **22**(12), pp. 2071-2083.
#'
#' Sahai, H. and Khurshid, A. (1993) Confidence intervals for the ratio of
#' two Poisson means. *The Mathematical Scientist*, **18**, pp. 43-50.
#'
#' @seealso [poissonCI()], [poissonDiffCI()], [binomRatioCI()],
#' [poisson.test()]
#'
#' @examples
#'
#' # 15 events in 100 person-years compared with
#' # 6 events in 120 person-years
#' poissonRatioCI(15, 100, 6, 120)
#'
#' # the exact interval corresponds to the two-sample conditional test
#' poisson.test(c(15, 6), c(100, 120))$conf.int
#'
#' poissonRatioCI(15, 100, 6, 120, method = "mid-p")
#'
#' # zero counts are handled explicitly
#' poissonRatioCI(0, 100, 6, 120)
#'
#' @concept confidence-interval
#' @concept rate
#'
#' @export
poissonRatioCI <- function(x1, n1 = 1, x2, n2 = 1, conf.level = 0.95,
                           sides = c("two.sided", "left", "right"),
                           method = c("exact", "mid-p", "wald-log")) {

  sides <- match.arg(sides)
  method <- match.arg(method)

  .poissonTwoSampleCI(x1, n1, x2, n2, conf.level, sides, method,
                      fun = function(x1, n1, x2, n2, conf.level, sides, method) {

    a <- if(sides == "two.sided") (1 - conf.level) / 2 else 1 - conf.level
    est <- (x1 / n1) / (x2 / n2)
    m <- x1 + x2
    f <- n2 / n1                       # maps the binomial odds back to the rate ratio

    if(m == 0) {

      # no events at all: the ratio 0/0 is undefined, though the data still
      # bound it by the whole positive half line
      est <- NA_real_
      lci <- 0
      uci <- Inf

    } else if(method == "wald-log") {

      if(x1 == 0 || x2 == 0) {
        lci <- 0
        uci <- Inf
      } else {
        h <- stats::qnorm(1 - a) * sqrt(1 / x1 + 1 / x2)
        lci <- exp(log(est) - h)
        uci <- exp(log(est) + h)
      }

    } else {

      if(method == "exact") {
        pl <- if(x1 == 0) 0 else stats::qbeta(a, x1, m - x1 + 1)
        pu <- if(x1 == m) 1 else stats::qbeta(1 - a, x1 + 1, m - x1)
        lci <- f * pl / (1 - pl)
        uci <- if(pu == 1) Inf else f * pu / (1 - pu)
      } else {
        # mid-p: half the point mass at the observed count enters each tail;
        # the limits come back as log-odds
        el <- if(x1 == 0) -Inf else .midpSolve(function(p)
                0.5 * stats::dbinom(x1, m, p) +
                stats::pbinom(x1, m, p, lower.tail = FALSE) - a)
        eu <- if(x1 == m) Inf else .midpSolve(function(p)
                0.5 * stats::dbinom(x1, m, p) +
                stats::pbinom(x1 - 1, m, p) - a)
        lci <- f * exp(el)
        uci <- f * exp(eu)
      }
    }

    if(sides == "left")  uci <- Inf
    if(sides == "right") lci <- 0

    c(est, lci, uci)
  })
}
