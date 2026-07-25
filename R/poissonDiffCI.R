
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


#' Confidence Intervals for the Difference of Two Poisson Rates
#'
#' Computes a confidence interval for the difference \eqn{x_1/n_1 - x_2/n_2} of
#' two independent Poisson rates.
#'
#' @details
#' \code{"mover"} (default) applies the method of variance estimates recovery of
#' Zou and Donner (2008): exact Garwood limits are obtained for each rate
#' separately and recombined, which keeps the interval close to nominal coverage
#' even for small counts. \code{"wald"} uses the asymptotic normal
#' approximation and degenerates to zero width when both counts are zero.
#'
#' Arguments are recycled to a common length. \code{sides} names the side
#' carrying the finite limit: \code{"left"} yields \eqn{[lci, \infty)} and
#' \code{"right"} yields \eqn{(-\infty, uci]}, with the full error probability
#' placed on that one side.
#'
#' @param x1 number of events in the first sample
#' @param n1 time base for the first event count
#' @param x2 number of events in the second sample
#' @param n2 time base for the second event count
#' @param conf.level confidence level, defaults to 0.95
#' @param sides a character string specifying the side of the confidence
#'   interval, one of \code{"two.sided"} (default), \code{"left"} or
#'   \code{"right"}
#' @param method character string specifying the method, either \code{"mover"}
#'   (default) or \code{"wald"}
#'
#' @return If recycling yields a single case, a named numeric vector with
#'   elements \code{est}, \code{lci} and \code{uci}. Otherwise a data frame with
#'   one row per case, whose first three columns are \code{est}, \code{lci} and
#'   \code{uci} and whose remaining columns contain the recycled argument
#'   values.
#'
#' @references
#' Zou, G. Y., Donner, A. (2008) Construction of confidence limits about effect
#' measures: a general approach. \emph{Statistics in Medicine}, \bold{27}(10),
#' 1693-1702.
#'
#' @seealso \code{\link{poissonCI}}, \code{\link{poissonRatioCI}},
#'   \code{\link{binomDiffCI}}
#'
#' @family ci.proportion
#' @concept rate
#'
#' @examples
#' ## 15 events in 100 person-years against 6 events in 120 person-years
#' poissonDiffCI(15, 100, 6, 120)
#' ## [1] est 0.10000000 lci 0.01155263 uci 0.20241565
#'
#' poissonDiffCI(15, 100, 6, 120, method = "wald")
#'
#' ## recycling returns one row per case
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


#' Confidence Intervals for the Ratio of Two Poisson Rates
#'
#' Computes a confidence interval for the ratio \eqn{(x_1/n_1) / (x_2/n_2)} of
#' two independent Poisson rates.
#'
#' @details
#' \code{"exact"} (default) conditions on the total count: given
#' \eqn{x_1 + x_2}, the count \eqn{x_1} is binomial with success probability
#' \eqn{n_1\theta/(n_1\theta + n_2)}, so a Clopper-Pearson interval for that
#' probability maps directly to one for the rate ratio \eqn{\theta}. This is the
#' construction underlying \code{\link[stats]{poisson.test}} for two samples.
#' \code{"midp"} replaces the Clopper-Pearson limits by their mid-p
#' counterparts, which are shorter and closer to nominal coverage at the cost of
#' exactness. \code{"wald-log"} is the asymptotic interval symmetric on the log
#' scale and is undefined when either count is zero, in which case
#' \eqn{(0, \infty)} is returned.
#'
#' Arguments are recycled to a common length. \code{sides} names the side
#' carrying the finite limit, and since the ratio is bounded below by zero the
#' open end of a one-sided interval is reported as \code{0} or \code{Inf}
#' rather than as \eqn{\pm\infty}.
#'
#' @param x1 number of events in the first sample
#' @param n1 time base for the first event count
#' @param x2 number of events in the second sample
#' @param n2 time base for the second event count
#' @param conf.level confidence level, defaults to 0.95
#' @param sides a character string specifying the side of the confidence
#'   interval, one of \code{"two.sided"} (default), \code{"left"} or
#'   \code{"right"}
#' @param method character string specifying the method, one of \code{"exact"}
#'   (default), \code{"midp"} or \code{"wald-log"}
#'
#' @return If recycling yields a single case, a named numeric vector with
#'   elements \code{est}, \code{lci} and \code{uci}. Otherwise a data frame with
#'   one row per case, whose first three columns are \code{est}, \code{lci} and
#'   \code{uci} and whose remaining columns contain the recycled argument
#'   values.
#'
#' @references
#' Sahai, H., Khurshid, A. (1993) Confidence intervals for the ratio of two
#' Poisson means. \emph{The Mathematical Scientist}, \bold{18}, 43-50.
#'
#' Graham, P. L., Mengersen, K., Morton, A. P. (2003) Confidence limits for the
#' ratio of two rates based on likelihood scores. \emph{Statistics in Medicine},
#' \bold{22}(12), 2071-2083.
#'
#' @seealso \code{\link{poissonCI}}, \code{\link{poissonDiffCI}},
#'   \code{\link{binomRatioCI}}, \code{\link[stats]{poisson.test}}
#'
#' @family ci.proportion
#' @concept rate
#'
#' @examples
#' ## 15 events in 100 person-years against 6 events in 120 person-years
#' poissonRatioCI(15, 100, 6, 120)
#' ## [1] est 3.000000 lci 1.099947 uci 9.437411
#'
#' ## agrees with the two-sample conditional test in stats
#' poisson.test(c(15, 6), c(100, 120))$conf.int
#'
#' poissonRatioCI(15, 100, 6, 120, method = "midp")
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
