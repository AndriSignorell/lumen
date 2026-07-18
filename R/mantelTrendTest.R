
#' Mantel Linear-by-Linear Association Test
#'
#' A chi-squared test for linear association between two ordinal variables
#' in a two-way contingency table, using row and column scores.
#'
#' Tests for linear trend in a two-way \eqn{r \times c}{r x c} contingency
#' table by computing the Mantel linear-by-linear association statistic
#' \deqn{Q_{MH} = (n - 1) \cdot r^2}{Q = (n-1) * r^2}
#' where \eqn{r} is the Pearson correlation between the row and column
#' variables using the supplied scores. Under the null hypothesis of no
#' linear association, \eqn{Q} has an asymptotic chi-squared distribution
#' with one degree of freedom.
#'
#' This test is sometimes called the Mantel-Haenszel chi-squared test for
#' trend, but it is \emph{not} the stratified Cochran-Mantel-Haenszel test
#' for \eqn{2 \times 2 \times k}{2 x 2 x k} tables (see
#' \code{\link[stats]{mantelhaen.test}} for that). It is a score test for
#' ordinal association, also known as the linear-by-linear association
#' test.
#'
#' Both variables should be measured on an ordinal scale. If \code{x} has
#' numeric row and/or column \code{dimnames} (e.g. income brackets or
#' dose levels stored as label strings that parse as numbers), those
#' values are used as the default scores; otherwise the default is
#' \code{1:nrow(x)} resp. \code{1:ncol(x)}, i.e. the categories are
#' assumed equally spaced. This mirrors the scoring convention used by
#' \code{\link{cochranArmitageTest}}. The choice of scores affects the
#' result: any monotone scores are permitted; non-monotone scores
#' (neither strictly increasing nor strictly decreasing) produce a
#' warning, since \eqn{r} would then no longer reflect a consistent
#' ordinal trend.
#'
#' @param x a numeric matrix of counts (\eqn{r \times c}{r x c})
#' @param srow numeric vector of row scores; length must equal
#' \code{nrow(x)}. Defaults to the numeric row \code{dimnames} of
#' \code{x} if present, otherwise \code{1:nrow(x)}. See the Details.
#' @param scol numeric vector of column scores; length must equal
#' \code{ncol(x)}. Defaults to the numeric column \code{dimnames} of
#' \code{x} if present, otherwise \code{1:ncol(x)}. See the Details.
#'
#' @return A list of class \code{"htest"} containing:
#' \item{statistic}{the Mantel linear association chi-squared statistic}
#' \item{parameter}{degrees of freedom (always 1)}
#' \item{p.value}{the p-value}
#' \item{estimate}{the Pearson correlation coefficient \eqn{r}}
#' \item{method}{a character string describing the test}
#' \item{data.name}{a character string giving the name of the data}
#'
#' @seealso \code{\link[stats]{mantelhaen.test}} for the stratified
#' Cochran-Mantel-Haenszel test, \code{\link{chisq.test}} for the general
#' chi-squared test of independence, \code{\link{cochranArmitageTest}} for
#' a related trend test with a binary response
#'
#' @references
#' Agresti, A. (2002) \emph{Categorical Data Analysis}, John Wiley & Sons,
#' pp. 57, 86.
#'
#' Mantel, N. (1963) Chi-square tests with one degree of freedom:
#' extensions of the Mantel-Haenszel procedure. \emph{Journal of the
#' American Statistical Association}, 58, 690-700.
#'
#' @examples
#' ## Agresti (2002, p. 57) Job Satisfaction
#' Job <- matrix(c(1,2,1,0, 3,3,6,1, 10,10,14,9, 6,7,12,11), 4, 4,
#'               dimnames = list(
#'                 income       = c("< 15k","15-25k","25-40k","> 40k"),
#'                 satisfaction = c("VeryD","LittleD","ModerateS","VeryS")))
#'
#' mantelTrendTest(Job)
#' mantelTrendTest(Job, srow = c(7.5, 20, 32.5, 60))
#'
#' ## Automatic scores from numeric dimnames
#' dose <- matrix(c(10, 9, 10, 7, 0, 1, 0, 3), nrow = 4,
#'                dimnames = list(dose = c("0", "1", "2", "3"),
#'                                resp = c("no", "yes")))
#' mantelTrendTest(dose)  # srow taken as c(0, 1, 2, 3), not 1:4
#'
#' @family test.categorical
#' @concept categorical-test
#' @concept trend-test
#'
#' @export
mantelTrendTest <- function(x, srow = scores(x, MARGIN = 1L, method = "table"),
                            scol = scores(x, MARGIN = 2L, method = "table")) {

  DNAME <- deparse1(substitute(x))

  if (!is.numeric(x) || length(dim(x)) != 2L)
    stop("'x' must be a numeric matrix")

  if (any(x < 0, na.rm = TRUE) || any(!is.finite(x)))
    stop("all entries of 'x' must be nonnegative and finite")

  if (sum(x) <= 1)
    stop("'x' must contain at least 2 observations")

  if (length(srow) != nrow(x))
    stop("'srow' must have the same length as nrow(x)")

  if (length(scol) != ncol(x))
    stop("'scol' must have the same length as ncol(x)")

  if (!(all(diff(srow) > 0) || all(diff(srow) < 0)))
    warning("'srow' is neither strictly increasing nor strictly ",
            "decreasing; scores should be monotone and ordinal",
            call. = FALSE)

  if (!(all(diff(scol) > 0) || all(diff(scol) < 0)))
    warning("'scol' is neither strictly increasing nor strictly ",
            "decreasing; scores should be monotone and ordinal",
            call. = FALSE)

  r <- .pearsonCor(x, srow = srow, scol = scol)
  STATISTIC <- (sum(x) - 1) * r^2

  structure(
    list(
      statistic = c("X-squared" = STATISTIC),
      parameter = c(df = 1L),
      p.value   = pchisq(STATISTIC, df = 1L, lower.tail = FALSE),
      estimate  = c(r = r),
      method    = "Mantel linear-by-linear association test",
      data.name = DNAME
    ),
    class = "htest"
  )
}



# == internal helper functions ============================================

.pearsonCor <- function(x, srow = 1:nrow(x), scol = 1:ncol(x)) {

  n    <- sum(x)
  ubar <- sum(rowSums(x) * srow) / n
  vbar <- sum(colSums(x) * scol) / n

  ssr <- sum(rowSums(x) * (srow - ubar)^2)
  ssc <- sum(colSums(x) * (scol - vbar)^2)

  if (ssr <= 0 || ssc <= 0)
    stop("row or column scores have zero variance; ",
         "check that scores are not all identical and that ",
         "more than one row/column is occupied", call. = FALSE)

  ssrc <- sum(x * outer(srow - ubar, scol - vbar))

  ssrc / sqrt(ssr * ssc)
}
