
#' Games-Howell Post-Hoc Test for Pairwise Mean Comparisons Under Unequal Variances
#'
#' Performs all pairwise comparisons between group means without assuming equal
#' variances, using a separate Welch-type standard error and Satterthwaite
#' degrees of freedom for every pair and referring the result to the
#' studentized range distribution.
#'
#' @details
#' Tukey's HSD, \code{\link{scheffeTest}} and the parametric methods in
#' \code{\link{postHoc}} all rely on a single pooled error variance. When the
#' group variances differ, that pooled estimate is wrong for every pair that
#' does not happen to match it, and the procedure loses its nominal level -
#' liberally when the larger variance sits in the smaller group, conservatively
#' in the opposite case. Games and Howell (1976) replace the pooled term by the
#' pairwise Welch standard error, which makes this the post-hoc counterpart of
#' the Welch t test and of \code{\link[stats]{oneway.test}}.
#' \code{\link{yuenTTest}} addresses heteroscedasticity and non-normality with
#' trimmed means and separate Winsorized variance estimates, but provides no
#' all-pairs multiple-comparison procedure, so it is not the two-sample form of
#' this procedure.
#'
#' For groups \eqn{i} and \eqn{j} the statistic is
#' \deqn{q = \frac{|\bar{x}_j - \bar{x}_i|}{\sqrt{(s_i^2/n_i + s_j^2/n_j)/2}},}
#' referred to the studentized range distribution with \eqn{k} means and
#' Satterthwaite degrees of freedom. Because \eqn{q} already accounts for the
#' number of groups, the p-values and interval widths are simultaneous over all
#' \eqn{k(k-1)/2} comparisons; no further multiplicity adjustment is applied or
#' needed.
#'
#' With equal group sizes and equal variances the standard error reduces exactly
#' to Tukey's, and the two procedures differ only in the degrees of freedom -
#' pairwise \eqn{2(n-1)} instead of pooled \eqn{k(n-1)} - which makes
#' Games-Howell slightly the more conservative of the two in that case.
#'
#' Every group needs at least two non-missing observations, since each
#' contributes its own variance estimate. Two obstacles can still leave an
#' individual pair incomputable: both groups constant, which leaves no scale to
#' studentize by, and Welch degrees of freedom below two, where the studentized
#' range is undefined. The latter is not an exotic case - with only two
#' observations in each of two groups the Welch degrees of freedom lie in
#' \eqn{[1, 2]} and reach two only when the two variances are exactly equal.
#' Affected pairs are reported as \code{NA} with a warning while the remaining
#' comparisons are kept. Missing values are removed casewise.
#'
#' @param x a numeric vector of observations, an \code{aov} object, or a formula
#'   of the form \code{lhs ~ rhs} with a numeric \code{lhs} and a grouping
#'   \code{rhs}
#' @param g a vector or factor giving the group for each element of \code{x}
#' @param conf.level confidence level of the simultaneous intervals
#' @param formula a formula of the form \code{lhs ~ rhs}
#' @param data an optional data frame containing the model variables
#' @param subset an optional vector specifying a subset of observations
#' @param na.action a function indicating what should happen when the data
#'   contain \code{NA}s
#' @param \dots further arguments, passed to the default method
#'
#' @return An object of class \code{"PostHocTest"}: a list with one matrix,
#'   named after the grouping variable. The matrix has columns \code{diff} for
#'   the observed mean difference (second group minus first), \code{lci} and
#'   \code{uci} for the simultaneous confidence limits, and \code{pval} for the
#'   simultaneous p-value. Print and plot methods are available for class
#'   \code{"PostHocTest"}.
#'
#' @references
#' Games, P. A., Howell, J. F. (1976) Pairwise multiple comparison procedures
#' with unequal n's and/or variances: a Monte Carlo study. \emph{Journal of
#' Educational Statistics}, \bold{1}(2), 113-125.
#'
#' @seealso [TukeyHSD()], [yuenTTest]
#'
#' @family test.posthoc
#' @concept multiple-comparison
#'
#' @examples
#' gamesHowellTest(breaks ~ tension, data = warpbreaks)
#'
#' ## compare with the pooled-variance procedures
#' gamesHowellTest(breaks ~ tension, data = warpbreaks)$tension[, "pval"]
#' ## [1] 0.080082272 0.006355163 0.251298244
#'
#' TukeyHSD(aov(breaks ~ tension, data = warpbreaks))
#'
#' @export
gamesHowellTest <- function(x, ...) UseMethod("gamesHowellTest")


#' @rdname gamesHowellTest
#' @export
gamesHowellTest.default <- function(x, g, conf.level = 0.95, ...) {

  gname <- deparse1(substitute(g))

  if(!is.numeric(x))
    stop("'x' must be a numeric vector")
  if(length(x) != length(g))
    stop("'x' and 'g' must have the same length")
  if(!is.numeric(conf.level) || length(conf.level) != 1L || is.na(conf.level) ||
     conf.level <= 0 || conf.level >= 1)
    stop("'conf.level' must be a single value in (0, 1)")

  ok <- stats::complete.cases(x, g)
  x <- as.numeric(x[ok])
  g <- droplevels(factor(g[ok]))

  # an infinite value makes var() NaN, which would turn 'degenerate' into NA
  # and make the if() below fail rather than warn
  if(any(!is.finite(x)))
    stop("'x' must contain only finite values")

  k <- nlevels(g)
  if(k < 2L)
    stop("'g' must have at least two levels with non-missing observations")

  n <- tapply(x, g, length)
  m <- tapply(x, g, mean)
  v <- tapply(x, g, stats::var)

  if(any(n < 2L))
    stop("every group needs at least 2 non-missing observations, ",
         "since each contributes its own variance estimate")

  pairs <- utils::combn(k, 2L)
  i <- pairs[1L, ]
  j <- pairs[2L, ]

  a <- v[i] / n[i]
  b <- v[j] / n[j]
  diff <- m[j] - m[i]

  # Tukey's studentized range is tabulated for the difference divided by a
  # standard error that carries the factor 1/2, hence (a + b) / 2 here.
  se <- sqrt((a + b) / 2)
  df <- (a + b)^2 / (a^2 / (n[i] - 1) + b^2 / (n[j] - 1))

  # Two obstacles make a pair incomputable rather than merely imprecise, and
  # both yield NA so that the remaining comparisons survive:
  # (1) both groups constant, leaving no scale to studentize by;
  # (2) fewer than 2 Welch degrees of freedom, where ptukey and qtukey are
  #     undefined and return NaN with a domain warning.
  degenerate <- !(a + b > 0)
  lowDf <- !degenerate & !(df >= 2)

  if(any(degenerate))
    warning(sum(degenerate), " of ", length(degenerate), " comparisons involve ",
            "two groups that are both constant, for which the statistic is ",
            "undefined; these are reported as NA", call. = FALSE)
  if(any(lowDf))
    warning(sum(lowDf), " of ", length(lowDf), " comparisons have fewer than 2 ",
            "Welch degrees of freedom, where the studentized range is not ",
            "defined; these are reported as NA", call. = FALSE)

  pval <- rep(NA_real_, length(diff))
  half <- rep(NA_real_, length(diff))
  ok <- !degenerate & !lowDf
  pval[ok] <- stats::ptukey(abs(diff[ok]) / se[ok], nmeans = k, df = df[ok],
                            lower.tail = FALSE)
  half[ok] <- stats::qtukey(conf.level, nmeans = k, df = df[ok]) * se[ok]

  out <- cbind(diff = diff, lci = diff - half, uci = diff + half, pval = pval)
  rownames(out) <- paste(levels(g)[j], levels(g)[i], sep = "-")

  res <- list(out)
  names(res) <- gname

  attr(res, "orig.call") <- sys.call()
  attr(res, "conf.level") <- conf.level
  attr(res, "ordered") <- FALSE
  attr(res, "method") <- "Games-Howell"
  class(res) <- "PostHocTest"
  res
}


#' @rdname gamesHowellTest
#' @export
gamesHowellTest.formula <- function(formula, data, subset, na.action = na.pass, ...) {

  subsetExpr <- if(!missing(subset)) substitute(subset) else NULL

  mf <- resolveFormula(formula = formula, data = data, subset = subsetExpr,
                       na.action = na.action,
                       allowed = c("two-sample-independent",
                                   "n-sample-independent"))

  res <- gamesHowellTest.default(x = mf$x, g = mf$group, ...)
  names(res) <- all.vars(formula)[2L]
  attr(res, "orig.call") <- sys.call()
  res
}


#' @rdname gamesHowellTest
#' @export
gamesHowellTest.aov <- function(x, conf.level = 0.95, ...) {
  
  .stopIfCovariates(x)
  # check: no Welch-type procedure is available for ancova
  
  # Counting factors would silently accept an ANCOVA and drop the covariate, so
  # the model terms themselves are checked instead.
  tl <- attr(stats::terms(x), "term.labels")
  if(length(tl) != 1L)
    stop("the Games-Howell procedure needs a one-way model, ",
         "without covariates or further terms")
  if(!is.null(stats::weights(x)))
    stop("weighted models are not supported")

  mf <- stats::model.frame(x)
  # an offset need not appear among the term labels and would be dropped silently
  if(!is.null(stats::model.offset(mf)))
    stop("models with an offset are not supported")

  g <- mf[[tl]]
  if(!is.factor(g) && !is.character(g))
    stop("the predictor must be a factor")

  res <- gamesHowellTest.default(x = stats::model.response(mf), g = g,
                                 conf.level = conf.level, ...)
  names(res) <- tl
  attr(res, "orig.call") <- sys.call()
  res
}
