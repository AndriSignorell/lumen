#' Confidence Intervals for the Difference of Two Paired Binomial Proportions
#'
#' Computes confidence intervals for the difference between two binomial
#' proportions measured on the same subjects (matched pairs), e.g. the
#' accuracies of two classifiers on the same observations, or a binary outcome
#' before and after a treatment.
#'
#' The data form a paired 2x2 table
#' \tabular{lcc}{
#'                  \tab second: success \tab second: failure \cr
#'   first: success \tab \eqn{n_{11}}    \tab \eqn{n_{12}}    \cr
#'   first: failure \tab \eqn{n_{21}}    \tab \eqn{n_{22}}
#' }
#' with \eqn{N} pairs. The difference of the two success proportions is
#' estimated by
#' \deqn{
#' \hat{\delta} = \hat{p}_{1+} - \hat{p}_{+1} = \frac{n_{12} - n_{21}}{N},
#' }
#' so only the discordant pairs \eqn{n_{12}} and \eqn{n_{21}} move the
#' estimate, while the concordant pairs determine its precision. Treating the
#' two proportions as independent ([binomDiffCI()]) ignores this correlation
#' and gives intervals of the wrong width, usually far too wide.
#'
#' **Tango asymptotic score** (default):
#' Inverts the score test for \eqn{\delta}, with the nuisance parameter
#' estimated by restricted maximum likelihood under each tested value
#' (Tango, 1998). The paired counterpart of the Miettinen-Nurminen interval.
#'
#' **Newcombe hybrid score**:
#' Combines the Wilson score intervals of the two marginal proportions,
#' accounting for their correlation (Newcombe, 1998, method 10). The
#' correlation estimate is corrected towards zero for small positive values,
#' as proposed there.
#'
#' **Bonett-Price**:
#' A Wald interval centred at the Laplace-adjusted estimates of the
#' discordant probabilities, \eqn{(n_{12} + 1)/(N + 2)} and
#' \eqn{(n_{21} + 1)/(N + 2)} (Bonett and Price, 2012).
#'
#' **Agresti-Min**:
#' A Wald interval after adding 1/2 to each of the four cells (Agresti and
#' Min, 2005).
#'
#' **Wald** and **Wald with continuity correction**:
#' The large-sample normal approximation, the corrected version with the
#' discordance reduced by one in the variance, as in Fagerland et al.
#' (2014). Both are degenerate (zero width) without discordant pairs and are
#' not recommended.
#'
#' Limits outside the admissible range \eqn{[-1, 1]} are truncated.
#'
#' **Which interval should be used?**\cr
#' Fagerland, Lydersen and Laake (2014) recommend the Tango asymptotic score,
#' the Newcombe hybrid score and the Bonett-Price interval for general use;
#' the Wald intervals perform poorly.
#'
#' @param x a 2x2 table or matrix of paired counts as above (the first
#'   measurement in the rows, the second in the columns, success first), or a
#'   logical or 0/1 vector of the first measurement's outcomes when `y` is
#'   given.
#' @param y `NULL` (default) if `x` is a table, else a logical or 0/1 vector
#'   of the second measurement's outcomes, of the same length as `x` and with
#'   the pairs in the same order.
#' @param conf.level confidence level, default is 0.95.
#' @param sides a character string specifying the type of confidence interval:
#'   `"two.sided"` (default), `"left"`, or `"right"`.
#'   Partial matching is allowed.
#' @param method one or several of `"tango"`, `"newcombe-score"`,
#'   `"bonett-price"`, `"agresti-min"`, `"wald"`, `"wald-cc"`.
#' @param na.rm logical, whether pairs with a missing value in `x` or `y` are
#'   removed. With the default `FALSE` such pairs make the result `NA`.
#'
#' @return If a single method and level are requested, a named numeric vector
#' with elements:
#' \describe{
#'   \item{`est`}{point estimate of the difference of the paired proportions,
#'     \eqn{(n_{12} - n_{21})/N}.}
#'   \item{`lci`}{lower confidence interval bound.}
#'   \item{`uci`}{upper confidence interval bound.}
#' }
#'
#' Otherwise a data frame with one row per method and level; its first three
#' columns are `est`, `lci`, and `uci`, the remaining ones contain the
#' recycled argument values.
#'
#' @references
#' Agresti A, Min Y (2005).
#' Simple improved confidence intervals for comparing matched proportions.
#' *Statistics in Medicine*, 24, 729-740.
#'
#' Bonett DG, Price RM (2012).
#' Adjusted Wald confidence interval for a difference of binomial proportions
#' based on paired data.
#' *Journal of Educational and Behavioral Statistics*, 37, 479-488.
#'
#' Fagerland MW, Lydersen S, Laake P (2014).
#' Recommended tests and confidence intervals for paired binomial proportions.
#' *Statistics in Medicine*, 33, 2850-2875.
#'
#' Newcombe RG (1998).
#' Improved confidence intervals for the difference between binomial
#' proportions based on paired data.
#' *Statistics in Medicine*, 17, 2635-2650.
#'
#' Tango T (1998).
#' Equivalence test and confidence interval for the difference in proportions
#' for the paired-sample design.
#' *Statistics in Medicine*, 17, 891-908.
#'
#' @examples
#' # Airway hyper-responsiveness before and after stem cell transplantation
#' # (Bentur et al., 2009, from Fagerland et al., 2014)
#' tab <- matrix(c(1, 7, 1, 12), nrow = 2,
#'               dimnames = list(before = c("yes", "no"),
#'                               after  = c("yes", "no")))
#' binomPairedDiffCI(tab)
#'
#' binomPairedDiffCI(tab, method = c("tango", "newcombe-score",
#'                                   "bonett-price", "agresti-min",
#'                                   "wald", "wald-cc"))
#'
#' # accuracies of two classifiers on the same observations
#' set.seed(1)
#' ok1 <- rbinom(200, 1, 0.85) == 1
#' ok2 <- ifelse(runif(200) < 0.8, ok1, rbinom(200, 1, 0.8) == 1)
#' binomPairedDiffCI(ok1, ok2)
#'
#' # one-sided: is the first at least as accurate as the second?
#' binomPairedDiffCI(ok1, ok2, conf.level = 0.95, sides = "left")
#'
#' @seealso [binomDiffCI()] for independent proportions,
#'   [mcnemar.test()] for the corresponding test
#' @family ci.proportion
#' @concept confidence-interval
#' @concept proportion
#' @concept paired-samples
#'
#' @export
binomPairedDiffCI <- function(x, y = NULL,
                              conf.level = 0.95,
                              sides = c("two.sided", "left", "right"),
                              method = c("tango",
                                         "newcombe-score",
                                         "bonett-price",
                                         "agresti-min",
                                         "wald",
                                         "wald-cc"),
                              na.rm = FALSE) {

  sides <- match.arg(sides)

  method <- if (missing(method))
    eval(formals(sys.function())$method)[1L]
  else
    .resolveMethod(method, several.ok = TRUE)

  checkFlag(na.rm)
  for (cl in conf.level) .checkSidedLevel(cl, sides)

  n <- .pairedTable(x, y, na.rm)

  res <- .recycleApply(.binomPairedDiffCI_engine,
                       n11 = n[1L, 1L], n12 = n[1L, 2L],
                       n21 = n[2L, 1L], n22 = n[2L, 2L],
                       conf.level = conf.level,
                       sides = sides,
                       method = method)

  if (length(res) == 1L)
    return(res[[1L]])

  rc <- as.data.frame(attr(res, "recycle"))
  data.frame(do.call(rbind, res),
             rc[, c("conf.level", "sides", "method")])
}



.binomPairedDiffCI_engine <- function(n11, n12, n21, n22,
                                      conf.level, sides, method) {

  N <- n11 + n12 + n21 + n22
  est <- (n12 - n21) / N

  # an undefined table (missing values not removed) gives an NA interval
  if (is.na(N))
    return(c(est = NA_real_, lci = NA_real_, uci = NA_real_))

  alpha <- .sidesAlpha(conf.level, sides)
  z <- stats::qnorm(1 - alpha / 2)

  CI <- switch(method,
               "tango"          = .bpdci.tango(n12, n21, N, z),
               "newcombe-score" = .bpdci.newcombe(n11, n12, n21, n22, alpha),
               "bonett-price"   = .bpdci.bonettPrice(n12, n21, N, z),
               "agresti-min"    = .bpdci.agrestiMin(n11, n12, n21, n22, z),
               "wald"           = .bpdci.wald(n12, n21, N, z, correct = FALSE),
               "wald-cc"        = .bpdci.wald(n12, n21, N, z, correct = TRUE),
               stop(gettextf("Unknown method '%s'.", method)))

  # clamping to [-1, 1] and opening the free side, design_rules 8.1.4
  c(est = est, applySides(unname(CI), sides, lo = -1, hi = 1))
}



# ===============================================================
# internal helper functions


# The paired 2x2 table from a table/matrix or from two outcome vectors.
# A missing value without na.rm gives a table of NA counts, which the
# engine turns into an NA interval (design_rules 3.4).
.pairedTable <- function(x, y, na.rm) {

  if (is.null(y)) {

    if (!(is.matrix(x) || is.table(x)) || !identical(dim(x), c(2L, 2L)))
      stop("'x' must be a 2x2 table or matrix if 'y' is not given",
           call. = FALSE)

    n <- matrix(as.numeric(x), 2L, 2L)
    if (anyNA(n))
      return(matrix(NA_real_, 2L, 2L))
    if (any(n < 0) || any(n != round(n)) || any(!is.finite(n)))
      stop("'x' must contain non-negative integer counts", call. = FALSE)

  } else {

    if (length(x) != length(y))
      stop("'x' and 'y' must have the same length", call. = FALSE)

    for (v in list(x, y))
      if (!(is.logical(v) || (is.numeric(v) && all(v %in% c(0, 1, NA)))))
        stop("'x' and 'y' must be logical or 0/1 vectors", call. = FALSE)

    # paired data: a pair with a missing value is incomplete as a whole
    ok <- !is.na(x) & !is.na(y)
    if (!all(ok) && !na.rm)
      return(matrix(NA_real_, 2L, 2L))

    x <- as.logical(x[ok])
    y <- as.logical(y[ok])
    n <- matrix(c(sum(x & y), sum(!x & y), sum(x & !y), sum(!x & !y)), 2L, 2L)
  }

  if (sum(n) == 0)
    stop("there are no pairs to compute an interval from", call. = FALSE)

  n
}


.bpdci.wald <- function(n12, n21, N, z, correct) {

  est <- (n12 - n21) / N
  d <- if (correct) max(abs(n12 - n21) - 1, 0) else n12 - n21
  SE <- sqrt(max((n12 + n21) - d^2 / N, 0)) / N

  c(lci = est - z * SE, uci = est + z * SE)
}


.bpdci.agrestiMin <- function(n11, n12, n21, n22, z) {

  m12 <- n12 + 0.5
  m21 <- n21 + 0.5
  M <- n11 + n12 + n21 + n22 + 2

  est <- (m12 - m21) / M
  SE <- sqrt((m12 + m21) - (m12 - m21)^2 / M) / M

  c(lci = est - z * SE, uci = est + z * SE)
}


.bpdci.bonettPrice <- function(n12, n21, N, z) {

  p12 <- (n12 + 1) / (N + 2)
  p21 <- (n21 + 1) / (N + 2)

  est <- p12 - p21
  SE <- sqrt((p12 + p21 - (p12 - p21)^2) / (N + 2))

  c(lci = est - z * SE, uci = est + z * SE)
}


.bpdci.newcombe <- function(n11, n12, n21, n22, alpha) {

  N <- n11 + n12 + n21 + n22
  n1p <- n11 + n12; n2p <- n21 + n22      # first measurement
  np1 <- n11 + n21; np2 <- n12 + n22      # second measurement

  p1 <- n1p / N
  p2 <- np1 / N
  est <- p1 - p2

  w1 <- binomCI(x = n1p, n = N, conf.level = 1 - alpha, method = "wilson")
  w2 <- binomCI(x = np1, n = N, conf.level = 1 - alpha, method = "wilson")
  l1 <- w1[["lci"]]; u1 <- w1[["uci"]]
  l2 <- w2[["lci"]]; u2 <- w2[["uci"]]

  # correlation of the two marginal proportions, with Newcombe's correction
  # of small positive values towards zero
  psi <- if (n1p == 0 || n2p == 0 || np1 == 0 || np2 == 0) {
    0
  } else {
    A <- n11 * n22 - n12 * n21
    nprod <- sqrt(n1p * n2p * np1 * np2)
    if (A > N / 2) (A - N / 2) / nprod
    else if (A >= 0) 0
    else A / nprod
  }

  c(lci = est - sqrt((p1 - l1)^2 + (u2 - p2)^2 - 2 * psi * (p1 - l1) * (u2 - p2)),
    uci = est + sqrt((p2 - l2)^2 + (u1 - p1)^2 - 2 * psi * (p2 - l2) * (u1 - p1)))
}


# Tango's score interval: the deltas whose score statistic stays within
# +-z. The statistic decreases in delta, so each bound is the single root of
# T(delta) -+ z on (-1, 1); a bound at the boundary of the range is taken
# there.
.bpdci.tango <- function(n12, n21, N, z) {

  score <- function(delta) {
    # restricted ML estimate of p21 under delta
    A <- 2 * N
    B <- -n12 - n21 + (2 * N - n12 + n21) * delta
    C <- -n21 * delta * (1 - delta)
    p21 <- (sqrt(B^2 - 4 * A * C) - B) / (2 * A)
    (n12 - n21 - N * delta) / sqrt(N * (2 * p21 + delta * (1 - delta)))
  }

  tol <- 1e-10
  est <- (n12 - n21) / N

  bound <- function(target, atEdge, edge) {
    if (atEdge)
      return(edge)
    f <- function(d) score(d) - target
    lo <- -1 + tol
    hi <- 1 - tol
    if (sign(f(lo)) == sign(f(hi)))
      return(edge)
    stats::uniroot(f, c(lo, hi), tol = tol)$root
  }

  c(lci = bound(z, est == -1, -1),
    uci = bound(-z, est == 1, 1))
}
