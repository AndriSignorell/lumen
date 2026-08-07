
# Exact (Garwood) limits for a single Poisson rate at one-sided level a.
# Identical to poissonCI(method = "exact"), inlined here so that the recycling
# loop does not pay for a full htest-style call per row.
.poissonExactLo <- function(x, n, a) if(x == 0) 0 else stats::qchisq(a, 2 * x) / (2 * n)
.poissonExactHi <- function(x, n, a) stats::qchisq(1 - a, 2 * (x + 1)) / (2 * n)


# Root of a mid-p equation on (0, 1). At extreme one-sided levels the equation
# can have no root at all, in which case the boundary nearest to solving it is
# the honest answer and maps to a limit of 0 or Inf.
.midpSolve <- function(h) {
  lo <- .Machine$double.eps
  hi <- 1 - .Machine$double.eps
  flo <- h(lo)
  fhi <- h(hi)
  if(!is.finite(flo) || !is.finite(fhi) || flo * fhi > 0)
    return(if(abs(flo) <= abs(fhi)) 0 else 1)
  # the default tolerance (~1.2e-4) is amplified by the odds transform below
  stats::uniroot(h, c(lo, hi), tol = 1e-10)$root
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


#' Confidence Interval for the Difference Between Two Poisson Rates
#'
#' Estimates the difference between two independent Poisson event rates and
#' calculates a confidence interval.
#'
#' @param x1 non-negative integer event count or vector of counts for the first
#'   sample
#' @param n1 positive exposure associated with `x1`, such as observation time,
#'   person-time, or population at risk; may be a vector and defaults to 1
#' @param x2 non-negative integer event count or vector of counts for the second
#'   sample
#' @param n2 positive exposure associated with `x2`; may be a vector and
#'   defaults to 1
#' @param conf.level numeric confidence level between 0 and 1; defaults to 0.95
#' @param sides type of confidence interval: `"two.sided"`, `"left"`, or
#'   `"right"`; may be abbreviated
#' @param method method used to calculate the confidence interval: `"mover"`
#'   or `"wald"`; may be abbreviated and defaults to `"mover"`
#'
#' @return If the arguments identify a single result, a named numeric vector
#'   with elements:
#'   \describe{
#'     \item{`est`}{estimated rate difference}
#'     \item{`lci`}{lower confidence bound}
#'     \item{`uci`}{upper confidence bound}
#'   }
#'   Otherwise, a `data.frame` containing these three columns followed by the
#'   recycled values of `x1`, `n1`, `x2`, `n2`, and `conf.level`.
#'
#' @details
#' The function assumes two independent counts
#' \deqn{X_i \sim \mathrm{Poisson}(n_i\lambda_i), \quad i = 1, 2.}
#' The parameter of interest is the rate difference
#' \eqn{\Delta = \lambda_1 - \lambda_2}, estimated by
#' \eqn{\hat{\Delta} = x_1/n_1 - x_2/n_2}.
#'
#' The available confidence-interval methods are:
#' \describe{
#'   \item{`"mover"`}{the method of variance estimates recovery (MOVER),
#'     which combines separate exact Garwood limits for the two rates}
#'   \item{`"wald"`}{the normal-approximation interval based on the standard
#'     error \eqn{\sqrt{x_1/n_1^2 + x_2/n_2^2}}}
#' }
#' The Wald interval has zero width when both counts are zero.
#'
#' For `sides = "left"`, the function returns a lower one-sided confidence
#' bound and sets `uci` to `Inf`. For `sides = "right"`, it returns an upper
#' one-sided confidence bound and sets `lci` to `-Inf`.
#'
#' The numeric arguments are recycled to a common length only when their
#' lengths are compatible whole multiples. Incompatible lengths produce an
#' error. `sides` and `method` must each identify a single choice.
#'
#' @references
#' Zou, G. Y. and Donner, A. (2008). Construction of confidence limits about
#' effect measures: a general approach. \emph{Statistics in Medicine},
#' \bold{27}(10), 1693--1702.
#'
#' @seealso [poissonCI()], [poissonRatioCI()], [binomDiffCI()]
#'
#' @concept confidence-interval
#' @concept rate
#'
#' @examples
#' # 15 events in 100 person-years compared with
#' # 6 events in 120 person-years
#' poissonDiffCI(15, 100, 6, 120)
#'
#' poissonDiffCI(15, 100, 6, 120, method = "wald")
#'
#' # A 95% lower confidence bound for the rate difference
#' poissonDiffCI(15, 100, 6, 120, sides = "left")
#'
#' # Recycling returns one row per comparison
#' poissonDiffCI(x1 = c(15, 20), n1 = 100, x2 = 6, n2 = 120)
#'
#' @export
poissonDiffCI <- function(x1, n1 = 1, x2, n2 = 1, conf.level = 0.95,
                          sides = c("two.sided", "left", "right"),
                          method = c("mover", "wald")) {

  sides <- match.arg(sides)
  method <- match.arg(method)

  .poissonTwoSampleCI(x1, n1, x2, n2, conf.level, sides, method,
                      fun = function(x1, n1, x2, n2, conf.level, sides, method) {

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


#' Confidence Interval for the Ratio of Two Poisson Rates
#'
#' Estimates the ratio of two independent Poisson event rates and calculates a
#' confidence interval.
#'
#' @param x1 non-negative integer event count or vector of counts for the first
#'   sample
#' @param n1 positive exposure associated with `x1`, such as observation time,
#'   person-time, or population at risk; may be a vector and defaults to 1
#' @param x2 non-negative integer event count or vector of counts for the second
#'   sample
#' @param n2 positive exposure associated with `x2`; may be a vector and
#'   defaults to 1
#' @param conf.level numeric confidence level between 0 and 1; defaults to 0.95
#' @param sides type of confidence interval: `"two.sided"`, `"left"`, or
#'   `"right"`; may be abbreviated
#' @param method method used to calculate the confidence interval: `"exact"`,
#'   `"midp"`, or `"wald-log"`; may be abbreviated and defaults to `"exact"`
#'
#' @return If the arguments identify a single result, a named numeric vector
#'   with elements:
#'   \describe{
#'     \item{`est`}{estimated rate ratio}
#'     \item{`lci`}{lower confidence bound}
#'     \item{`uci`}{upper confidence bound}
#'   }
#'   Otherwise, a `data.frame` containing these three columns followed by the
#'   recycled values of `x1`, `n1`, `x2`, `n2`, and `conf.level`.
#'
#' @details
#' The function assumes two independent counts
#' \deqn{X_i \sim \mathrm{Poisson}(n_i\lambda_i), \quad i = 1, 2.}
#' The parameter of interest is the rate ratio
#' \eqn{\theta = \lambda_1/\lambda_2}, estimated by
#' \eqn{\hat{\theta} = (x_1/n_1)/(x_2/n_2)}.
#'
#' The available confidence-interval methods are:
#' \describe{
#'   \item{`"exact"`}{the exact conditional interval obtained by conditioning
#'     on \eqn{x_1 + x_2} and transforming a Clopper--Pearson interval for the
#'     resulting binomial probability; this is the construction used by
#'     [stats::poisson.test()] for two samples}
#'   \item{`"midp"`}{the corresponding conditional mid-p interval, which is
#'     generally shorter but does not guarantee conservative coverage}
#'   \item{`"wald-log"`}{the asymptotic Wald interval, symmetric on the
#'     logarithmic scale}
#' }
#'
#' For `sides = "left"`, the function returns a lower one-sided confidence
#' bound and sets `uci` to `Inf`. For `sides = "right"`, it returns an upper
#' one-sided confidence bound and sets `lci` to 0, the lower limit of the
#' parameter space.
#'
#' The log-Wald interval cannot be calculated when either count is zero and is
#' then returned as `[0, Inf]`. If both counts are zero, the rate ratio and its
#' point estimate are undefined; the function returns `NA` for `est` and
#' `[0, Inf]` for the confidence interval for every method. If only `x2` is
#' zero, the point estimate is `Inf`.
#'
#' The numeric arguments are recycled to a common length only when their
#' lengths are compatible whole multiples. Incompatible lengths produce an
#' error. `sides` and `method` must each identify a single choice.
#'
#' @references
#' Sahai, H. and Khurshid, A. (1993). Confidence intervals for the ratio of two
#' Poisson means. \emph{The Mathematical Scientist}, \bold{18}, 43--50.
#'
#' Graham, P. L., Mengersen, K. and Morton, A. P. (2003). Confidence limits for
#' the ratio of two rates based on likelihood scores. \emph{Statistics in
#' Medicine}, \bold{22}(12), 2071--2083.
#'
#' @seealso [poissonCI()], [poissonDiffCI()], [binomRatioCI()],
#'   [stats::poisson.test()]
#'
#' @concept confidence-interval
#' @concept rate
#'
#' @examples
#' # 15 events in 100 person-years compared with
#' # 6 events in 120 person-years
#' poissonRatioCI(15, 100, 6, 120)
#'
#' # The exact interval agrees with the two-sample conditional test in stats
#' poisson.test(c(15, 6), c(100, 120))$conf.int
#'
#' poissonRatioCI(15, 100, 6, 120, method = "midp")
#'
#' # Zero counts are handled explicitly
#' poissonRatioCI(0, 100, 6, 120)
#'
#' @export
poissonRatioCI <- function(x1, n1 = 1, x2, n2 = 1, conf.level = 0.95,
                           sides = c("two.sided", "left", "right"),
                           method = c("exact", "midp", "wald-log")) {

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
      } else {
        # mid-p: half the point mass at the observed count enters each tail
        pl <- if(x1 == 0) 0 else .midpSolve(function(p)
                0.5 * stats::dbinom(x1, m, p) +
                stats::pbinom(x1, m, p, lower.tail = FALSE) - a)
        pu <- if(x1 == m) 1 else .midpSolve(function(p)
                0.5 * stats::dbinom(x1, m, p) +
                stats::pbinom(x1 - 1, m, p) - a)
      }
      lci <- f * pl / (1 - pl)
      uci <- if(pu == 1) Inf else f * pu / (1 - pu)
    }

    if(sides == "left")  uci <- Inf
    if(sides == "right") lci <- 0

    c(est, lci, uci)
  })
}
