
#' Mood's Median Test for Comparing Medians Across Independent Groups
#'
#' Tests whether two or more independent samples come from populations with the
#' same median. All observations are dichotomised at the pooled median and the
#' resulting 2 x k table is tested for independence.
#'
#' @details
#' Not to be confused with \code{\link[stats]{mood.test}}, which is Mood's
#' two-sample test for a difference in \emph{scale}. The median test described
#' here has no base R implementation.
#'
#' The procedure is the k-sample counterpart of \code{\link{signTest}} and
#' shares its robustness and its modest power: only the position of each
#' observation relative to the pooled median enters the statistic, so all
#' information about distance is discarded. Its asymptotic relative efficiency
#' against the F test under normality is \eqn{2/\pi}.
#'
#' That low efficiency is the price of asking a narrow question, and the
#' alternatives ask different ones rather than the same one better:
#' \code{\link{brunnerMunzelTest}} tests the relative effect
#' \eqn{P(X < Y) + \frac{1}{2}P(X = Y) = \frac{1}{2}} and
#' \code{\link{vanWaerdenTest}} tests equality of the distributions against
#' normal-score location alternatives. Neither is a test of equal medians, so
#' they are not drop-in replacements. Use the median test when the median is
#' genuinely the quantity of interest, or when only the side of a threshold is
#' trustworthy, as with coarsely recorded or thresholded data.
#'
#' \strong{Observations equal to the pooled median.} With an even total sample
#' size the pooled median usually falls between two observations and the
#' question does not arise. Otherwise \code{ties} decides: \code{"below"}
#' (default) counts them with the lower group, following Conover, \code{"above"}
#' with the upper group, and \code{"drop"} removes them, which retains a
#' symmetric above-versus-below classification at the cost of a smaller
#' effective sample size. Dropping them does not systematically enlarge the
#' p-value; it can move it either way.
#' \code{"below"} and \code{"above"} correspond to \code{mid.score} values of
#' \code{"0"} and \code{"1"} in \code{coin::median_test}; \code{"drop"} has no
#' counterpart there, since \code{mid.score = "0.5"} scores the median
#' observations rather than removing them.
#'
#' \code{method = "exact"} replaces the chi-squared approximation by Fisher's
#' exact test on the same table and is advisable when expected counts are small;
#' no test statistic is reported in that case. The continuity correction applies
#' only to a 2 x 2 table, as in \code{\link[stats]{chisq.test}}.
#'
#' Missing values are removed casewise.
#'
#' @param x a numeric vector of observations, or a formula of the form
#'   \code{lhs ~ rhs} with a numeric \code{lhs} and a grouping \code{rhs}
#' @param g a vector or factor giving the group for each element of \code{x}
#' @param ties how to treat observations exactly equal to the pooled median, one
#'   of \code{"below"} (default), \code{"above"} or \code{"drop"}
#' @param method the test applied to the 2 x k table, either \code{"chisq"}
#'   (default) or \code{"exact"}
#' @param correct logical, whether to apply the continuity correction; ignored
#'   unless the table is 2 x 2
#' @param formula a formula of the form \code{lhs ~ rhs}
#' @param data an optional data frame containing the model variables
#' @param subset an optional vector specifying a subset of observations
#' @param na.action a function indicating what should happen when the data
#'   contain \code{NA}s
#' @param \dots further arguments, passed to the default method
#'
#' @return An object of class \code{"htest"} with components
#'   \item{statistic}{Pearson's chi-squared statistic, \code{NULL} for
#'     \code{method = "exact"}}
#'   \item{parameter}{the degrees of freedom, \code{NULL} for
#'     \code{method = "exact"}}
#'   \item{p.value}{the p-value}
#'   \item{estimate}{the group medians}
#'   \item{observed}{the 2 x k table of counts above and below the pooled median}
#'   \item{expected}{the expected counts under independence}
#'   \item{grand.median}{the pooled median at which the data were dichotomised}
#'   \item{method}{a character string describing the test}
#'   \item{data.name}{a character string giving the names of the data}
#'
#' @references
#' Mood, A. M. (1950) \emph{Introduction to the Theory of Statistics}.
#' McGraw-Hill, New York, pp. 394-399.
#'
#' Conover, W. J. (1999) \emph{Practical Nonparametric Statistics}, 3rd edition.
#' Wiley, New York, pp. 218-223.
#'
#' @seealso \code{\link{signTest}}, \code{\link{brunnerMunzelTest}},
#'   \code{\link{vanWaerdenTest}}, \code{\link[stats]{mood.test}}
#'
#' @family test.location
#' @concept median
#'
#' @examples
#' moodMedianTest(breaks ~ tension, data = warpbreaks)
#'
#' ## the pooled median at which the data were split, and the resulting table
#' r <- moodMedianTest(breaks ~ tension, data = warpbreaks)
#' r$grand.median
#' ## [1] 26
#' r$observed
#'
#' ## small tables: use the exact test on the same 2 x k table
#' moodMedianTest(breaks ~ tension, data = warpbreaks, method = "exact")$p.value
#' ## [1] 0.03153902
#'
#' @export
moodMedianTest <- function(x, ...) UseMethod("moodMedianTest")


#' @rdname moodMedianTest
#' @export
moodMedianTest.default <- function(x, g, ties = c("below", "above", "drop"),
                                   method = c("chisq", "exact"),
                                   correct = TRUE, ...) {

  ties <- match.arg(ties)
  method <- match.arg(method)

  dname <- paste(deparse1(substitute(x)), "and", deparse1(substitute(g)))

  if(!is.numeric(x))
    stop("'x' must be a numeric vector")
  if(length(x) != length(g))
    stop("'x' and 'g' must have the same length")
  if(!is.logical(correct) || length(correct) != 1L || is.na(correct))
    stop("'correct' must be a single non-missing logical value")

  ok <- stats::complete.cases(x, g)
  x <- as.numeric(x[ok])
  g <- droplevels(factor(g[ok]))

  k <- nlevels(g)
  if(k < 2L)
    stop("'g' must have at least two levels with non-missing observations")

  grandMedian <- stats::median(x)

  # "drop" removes the observations sitting exactly on the pooled median; the
  # other two conventions assign them to the lower or the upper cell.
  keep <- if(ties == "drop") x != grandMedian else rep(TRUE, length(x))
  above <- if(ties == "above") x >= grandMedian else x > grandMedian

  side <- factor(ifelse(above[keep], "above", "below"),
                 levels = c("above", "below"))
  observed <- table(side, g[keep])
  names(dimnames(observed)) <- c("", "")

  n <- sum(observed)
  if(n < 2L)
    stop("no observations left after removing values equal to the pooled median")
  if(any(colSums(observed) == 0L))
    stop("at least one group has no observations left")
  if(any(rowSums(observed) == 0L))
    stop("all observations fall on the same side of the pooled median, ",
         "so the table carries no information")

  expected <- outer(rowSums(observed), colSums(observed)) / n

  if(method == "exact") {

    statistic <- NULL
    parameter <- NULL
    pval <- stats::fisher.test(observed, ...)$p.value
    mstr <- gettextf("Mood's median test (%d groups, Fisher's exact test)", k)

  } else {

    if(any(expected < 5))
      warning("chi-squared approximation may be incorrect", call. = FALSE)

    # Yates' correction is defined for a 2 x 2 table only, as in chisq.test()
    yates <- correct && all(dim(observed) == 2L)
    dev <- abs(observed - expected)
    if(yates)
      dev <- pmax(0, dev - 0.5)

    statistic <- c("X-squared" = sum(dev^2 / expected))
    parameter <- c(df = k - 1L)
    pval <- stats::pchisq(statistic, df = parameter, lower.tail = FALSE)
    mstr <- gettextf("Mood's median test (%d groups, chi-squared approximation%s)",
                     k, if(yates) " with continuity correction" else "")
  }

  res <- list(statistic = statistic,
              parameter = parameter,
              p.value = unname(pval),
              estimate = stats::setNames(as.numeric(tapply(x, g, stats::median)),
                                         levels(g)),
              observed = observed,
              expected = expected,
              grand.median = grandMedian,
              method = mstr,
              data.name = dname)

  class(res) <- "htest"
  res
}


#' @rdname moodMedianTest
#' @export
moodMedianTest.formula <- function(formula, data, subset, na.action = na.pass, ...) {

  subsetExpr <- if(!missing(subset)) substitute(subset) else NULL

  mf <- resolveFormula(formula = formula, data = data, subset = subsetExpr,
                       na.action = na.action,
                       allowed = c("two-sample-independent",
                                   "n-sample-independent"))

  res <- moodMedianTest.default(x = mf$x, g = mf$group, ...)
  res$data.name <- mf$data.name
  res
}
