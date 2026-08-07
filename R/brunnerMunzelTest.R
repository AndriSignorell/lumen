
# Splits above this count trigger the Monte-Carlo permutation path when
# 'exact' is left at its default. ~1e6 splits run in roughly 0.2 s.
.bmMaxExact <- 1e6

# Hard ceiling for exact = TRUE. At roughly 5e6 splits per second this is a few
# minutes of work; beyond it, enumeration is almost certainly a mistake.
.bmMaxForcedExact <- 1e9


#' Brunner-Munzel Test for Comparing Stochastic Dominance Between Two Independent Samples
#'
#' Tests the nonparametric Behrens-Fisher hypothesis that a randomly drawn
#' observation from one sample is as likely to be smaller as it is to be larger
#' than a randomly drawn observation from the other, without assuming equal
#' shapes or equal variances in the two groups.
#'
#' The estimated quantity is the relative effect
#' \deqn{p = P(X < Y) + \tfrac{1}{2} P(X = Y),}
#' the probability that a random draw from \code{y} exceeds a random draw from
#' \code{x}, ties counted half. It is identical to the Mann-Whitney statistic
#' scaled to \eqn{[0, 1]}, i.e. \eqn{U / (n_1 n_2)}, and is often reported as
#' the common language effect size.
#'
#' @details
#' The Wilcoxon-Mann-Whitney test is a valid test of \eqn{p = 1/2} only under
#' the hypothesis that both samples come from the same distribution. When the
#' two distributions differ in shape or spread, its null distribution is wrong
#' and the test does not keep its level. Brunner and Munzel (2000) studentize
#' the rank statistic with a separate variance estimate per group and refer it
#' to a t distribution with Satterthwaite degrees of freedom, which is the rank
#' analogue of the Welch correction. \code{\link{yuenTTest}} plays the same role
#' among the parametric location tests.
#'
#' \strong{Direction of \code{alternative}.} The alternative is stated in terms
#' of \eqn{p}, not in terms of \code{x} against \code{y}. \code{"greater"}
#' therefore means \eqn{p > p_0}, that is, \code{y} tends to produce the larger
#' values. This follows the published definition of the statistic and the
#' reported estimate, but it is the reverse of \code{\link[stats]{t.test}} and
#' \code{\link[stats]{wilcox.test}}, where \code{"greater"} refers to \code{x}.
#'
#' \strong{Choice of \code{method}.} The t approximation is liberal in small
#' samples; below roughly ten observations per group the studentized permutation
#' test of Neubert and Brunner (2007) should be preferred. Permuting the
#' studentized statistic rather than the raw rank statistic is what keeps the
#' permutation test asymptotically valid when the distributions differ, so this
#' is not an ordinary permutation test on ranks. It targets \eqn{H_0: p = 1/2}:
#' it is exact in finite samples under exchangeability of the group labels, and
#' remains asymptotically valid under the weaker nonparametric Behrens-Fisher
#' null \eqn{p = 1/2} alone. Note that \code{exact = TRUE} means exact
#' enumeration of the permutation distribution, which is not the same as a
#' finite-sample exact test under the general null. Only \eqn{p_0 = 1/2} is
#' implemented; other values of \code{p0} require an approximate method.
#'
#' The two-sided permutation p-value is the proportion of splits with
#' \eqn{|T^*| \ge |T|}, which is the convention used by the \pkg{brunnermunzel}
#' package implementation based on Neubert and Brunner (2007). It is exact under
#' exchangeability, but it does not in general equal twice the smaller one-sided
#' tail, because the permutation distribution need not be symmetric when the
#' group sizes differ or ties are present.
#'
#' With \code{exact = NULL} all \eqn{\binom{n_1 + n_2}{n_1}} splits are
#' enumerated when there are at most \code{1e6} of them, and \code{nPerm}
#' Monte-Carlo resamples are drawn otherwise. Monte-Carlo p-values use the
#' \eqn{(1 + k) / (1 + B)} correction and are therefore never zero.
#'
#' \strong{Confidence interval.} The interval is the studentized Wald interval
#' for \eqn{p} and is reported for every \code{method}, including the
#' permutation methods, which supply a p-value but no interval of their own. It
#' is clipped to \eqn{[0, 1]}, and for a one-sided \code{alternative} the open
#' end is reported at the range limit rather than as \eqn{\pm\infty}.
#'
#' \strong{Non-overlapping samples.} If every observation in one sample is
#' smaller than every observation in the other, both variance components are
#' zero and neither approximation is defined. The permutation test is used
#' instead, with a warning. For tie-free data small enough to enumerate, this
#' returns the smallest attainable two-sided p-value,
#' \eqn{2 / \binom{n_1 + n_2}{n_1}}; with ties the two mirrored fully separated
#' splits need not both exist, so the attainable minimum can be smaller. No Wald
#' interval exists in that case and \code{conf.int} is \code{NA}: a zero
#' standard error reflects an empty variance estimate, not certainty about
#' \eqn{p}.
#'
#' The formula method accepts \code{lhs ~ rhs} with exactly two groups on the
#' right-hand side; more than two levels are rejected by
#' [bedrock::resolveFormula].
#'
#' Missing values are removed.
#'
#' @param x a numeric vector or ordered factor
#' @param y a numeric vector or ordered factor
#' @param formula a formula of the form \code{lhs ~ rhs} where \code{lhs} is
#'   numeric and \code{rhs} a factor with two levels
#' @param p0 the value of \eqn{p} under the null hypothesis; must be \code{0.5}
#'   for \code{method = "permutation"}
#' @param alternative a character string specifying the alternative hypothesis
#'   for \eqn{p}, one of \code{"two.sided"} (default), \code{"less"} or
#'   \code{"greater"}
#' @param conf.level confidence level of the interval
#' @param method the inference method, one of \code{"t"} (default, Satterthwaite
#'   t approximation), \code{"permutation"} (studentized permutation test) or
#'   \code{"normal"} (asymptotic normal approximation)
#' @param exact logical, whether to enumerate all splits instead of sampling
#'   them; \code{NULL} (default) decides by the number of splits. Ignored unless
#'   \code{method = "permutation"}
#' @param nPerm number of Monte-Carlo resamples used when the permutation
#'   distribution is not enumerated
#' @param data an optional data frame containing the model variables
#' @param subset an optional vector specifying a subset of observations
#' @param na.action a function indicating what should happen when the data
#'   contain \code{NA}s
#' @param \dots further arguments, passed to the default method
#'
#' @return An object of class \code{"htest"} with components
#'   \item{statistic}{the studentized Brunner-Munzel statistic}
#'   \item{parameter}{the Satterthwaite degrees of freedom, \code{NULL} for
#'     \code{method = "permutation"}}
#'   \item{p.value}{the p-value}
#'   \item{conf.int}{confidence interval for \eqn{p}}
#'   \item{estimate}{the estimated relative effect \eqn{\hat{p}}}
#'   \item{null.value}{the value of \eqn{p} under the null hypothesis}
#'   \item{stderr}{the standard error of \eqn{\hat{p}}}
#'   \item{alternative}{the alternative hypothesis}
#'   \item{method}{a character string describing the test}
#'   \item{data.name}{a character string giving the names of the data}
#'
#' @references
#' Brunner, E., Munzel, U. (2000). The nonparametric Behrens-Fisher problem:
#' asymptotic theory and a small-sample approximation. \emph{Biometrical
#' Journal}, \bold{42}(1), 17-25.
#'
#' Neubert, K., Brunner, E. (2007). A studentized permutation test for the
#' non-parametric Behrens-Fisher problem. \emph{Computational Statistics and
#' Data Analysis}, \bold{51}(10), 5192-5204.
#'
#' @seealso [wilcox.test()]
#'
#' @family test.location
#' @concept rank-statistic
#' @concept effect-size
#'
#' @examples
#' ## Brunner & Munzel (2000), pain scores
#' x <- c(1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 2, 4, 1, 1)
#' y <- c(3, 3, 4, 3, 1, 2, 3, 1, 1, 5, 4)
#'
#' brunnerMunzelTest(x, y)
#'
#' ## the estimate is the common language effect size
#' brunnerMunzelTest(x, y)$estimate
#' ## [1] 0.7889610
#'
#' ## small groups: prefer the studentized permutation test, which is
#' ## enumerated exactly whenever the number of splits allows it
#' a <- c(1, 1, 2, 3, 3)
#' b <- c(4, 2, 5, 1, 6, 3, 7)
#' brunnerMunzelTest(a, b, method = "permutation")$p.value
#' ## [1] 0.09469697
#'
#' ## formula interface
#' d <- data.frame(score = c(x, y),
#'                 grp = rep(c("a", "b"), c(length(x), length(y))))
#' brunnerMunzelTest(score ~ grp, data = d)
#'
#' @export
brunnerMunzelTest <- function(x, ...) UseMethod("brunnerMunzelTest")


#' @rdname brunnerMunzelTest
#' @export
brunnerMunzelTest.default <- function(x, y, p0 = 0.5,
                                      alternative = c("two.sided", "less", "greater"),
                                      conf.level = 0.95,
                                      method = c("t", "permutation", "normal"),
                                      exact = NULL, nPerm = 10000L, ...) {

  alternative <- match.arg(alternative)
  method <- match.arg(method)

  dname <- paste(deparse1(substitute(x)), "and", deparse1(substitute(y)))

  # Ordered factors are legitimate input for a rank test, but only if both carry
  # the same level order: identical labels in a different order would otherwise
  # be silently mapped to opposing integer codes. Unordered factors and
  # character vectors must not be coerced at all.
  xOrdered <- is.ordered(x)
  yOrdered <- is.ordered(y)
  if(xOrdered || yOrdered) {
    if(!xOrdered || !yOrdered || !identical(levels(x), levels(y)))
      stop("'x' and 'y' must both be ordered factors with identical levels")
    x <- as.integer(x)
    y <- as.integer(y)
  }
  if(!is.numeric(x) || !is.numeric(y))
    stop("'x' and 'y' must be numeric vectors or ordered factors")

  x <- x[!is.na(x)]
  y <- y[!is.na(y)]
  n1 <- length(x)
  n2 <- length(y)
  n <- n1 + n2

  if(n1 < 2L || n2 < 2L)
    stop("both samples must contain at least 2 non-missing observations")
  if(!is.numeric(p0) || length(p0) != 1L || is.na(p0) || p0 <= 0 || p0 >= 1)
    stop("'p0' must be a single value in (0, 1)")
  if(!is.numeric(conf.level) || length(conf.level) != 1L || is.na(conf.level) ||
     conf.level <= 0 || conf.level >= 1)
    stop("'conf.level' must be a single value in (0, 1)")
  if(!is.numeric(nPerm) || length(nPerm) != 1L || !is.finite(nPerm) ||
     nPerm < 1 || nPerm != floor(nPerm) || nPerm > .Machine$integer.max)
    stop("'nPerm' must be a single positive integer no larger than ",
         ".Machine$integer.max")
  if(!is.null(exact) && (!is.logical(exact) || length(exact) != 1L || is.na(exact)))
    stop("'exact' must be NULL or a single non-missing logical value")

  nPerm <- as.integer(nPerm)

  r <- rank(c(x, y))
  core <- bm_core_cpp(r[seq_len(n1)], r[n1 + seq_len(n2)])
  phat <- core[["phat"]]
  se <- core[["se"]]
  df <- core[["df"]]

  # Both variance components vanish when the samples do not overlap at all, and
  # when every observation is tied. The studentized statistic then degenerates
  # and only the permutation test remains well defined.
  degenerate <- !isTRUE(se > 0)
  if(degenerate && method != "permutation") {
    if(p0 != 0.5)
      stop("the variance estimator is zero; no approximation is available for 'p0' other than 0.5")
    warning("the variance estimator is zero (the samples do not overlap, or all ",
            "observations are tied); using a permutation test instead of ",
            "the ", method, " approximation, and no Wald confidence interval ",
            "is available", call. = FALSE)
    method <- "permutation"
    # deliberately no 'exact <- TRUE': the usual auto-switch below picks
    # enumeration when it is affordable and sampling when it is not, so a
    # default call on large separated samples degrades instead of erroring
  }

  statistic <- if(!degenerate) {
    (phat - p0) / se
  } else if(phat == p0) {
    0                                   # every observation tied
  } else {
    sign(phat - p0) * Inf               # samples do not overlap
  }

  if(method == "permutation") {

    if(p0 != 0.5)
      stop("'method = \"permutation\"' is implemented only for 'p0 = 0.5'")

    nSplit <- choose(n, n1)
    if(is.null(exact))
      exact <- nSplit <= .bmMaxExact

    if(exact) {
      if(nSplit > .bmMaxForcedExact)
        stop(gettextf(paste("exact enumeration would require %.0f splits;",
                            "use 'exact = FALSE' to sample the permutation",
                            "distribution instead"), nSplit))
      if(nSplit > .bmMaxExact)
        message(gettextf("enumerating %.0f splits, this may take a while", nSplit))
      cnt <- bm_perm_exact_cpp(sort(r), n1, statistic)
      pval <- cnt[[alternative]] / cnt[["n"]]
      mstr <- gettextf("Brunner-Munzel test (exact studentized permutation, %.0f splits)",
                       cnt[["n"]])
    } else {
      cnt <- bm_perm_mc_cpp(sort(r), n1, statistic, nPerm)
      pval <- (1 + cnt[[alternative]]) / (1 + cnt[["n"]])
      mstr <- gettextf("Brunner-Munzel test (studentized permutation, %.0f resamples)",
                       cnt[["n"]])
    }

  } else {

    if(!is.null(exact))
      warning("'exact' is ignored unless 'method = \"permutation\"'", call. = FALSE)
    if(method == "normal")
      df <- Inf

    pval <- switch(alternative,
                   two.sided = 2 * stats::pt(-abs(statistic), df = df),
                   less      = stats::pt(statistic, df = df),
                   greater   = stats::pt(statistic, df = df, lower.tail = FALSE))
    mstr <- switch(method,
                   t      = "Brunner-Munzel test (t approximation)",
                   normal = "Brunner-Munzel test (normal approximation)")
  }

  # Studentized Wald interval, clipped to the parameter range. The open end of a
  # one-sided interval is reported at the range limit, not as +/-Inf.
  # A zero Wald standard error is an empty variance estimate, not certainty:
  # reporting c(phat, phat) would claim a degenerate interval at 0 or 1.
  ci <- if(degenerate) {
    c(NA_real_, NA_real_)
  } else switch(alternative,
                two.sided = phat + c(-1, 1) * stats::qt(1 - (1 - conf.level) / 2, df = df) * se,
                less      = c(0, phat + stats::qt(conf.level, df = df) * se),
                greater   = c(phat - stats::qt(conf.level, df = df) * se, 1))
  ci <- pmin(pmax(ci, 0), 1)
  attr(ci, "conf.level") <- conf.level

  estName <- "P(X < Y) + 0.5 * P(X = Y)"
  res <- list(statistic = c(W = statistic),
              parameter = if(method != "permutation") c(df = df),
              p.value = pval,
              conf.int = ci,
              estimate = stats::setNames(phat, estName),
              null.value = stats::setNames(p0, estName),
              stderr = se,
              alternative = alternative,
              method = mstr,
              data.name = dname)

  class(res) <- "htest"
  res
}


#' @rdname brunnerMunzelTest
#' @export
brunnerMunzelTest.formula <- function(formula, data, subset, na.action = na.pass, ...) {

  # 'subset' must be captured here and handed to resolveFormula() as an
  # unevaluated expression; the call has to be direct, since do.call() would
  # evaluate that expression in the wrong frame (hotellingsT2Test regression).
  subsetExpr <- if(!missing(subset)) substitute(subset) else NULL

  mf <- resolveFormula(formula = formula, data = data, subset = subsetExpr,
                       na.action = na.action,
                       allowed = "two-sample-independent")

  # x + group is the canonical access path; mf$y is convenience only
  grp <- split(mf$x, mf$group)

  res <- brunnerMunzelTest.default(x = grp[[1L]], y = grp[[2L]], ...)
  res$data.name <- mf$data.name
  res
}
