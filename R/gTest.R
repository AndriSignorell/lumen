
#' G-Test for Count Data
#'
#' A goodness-of-fit or test of independence based on the log-likelihood
#' ratio (G-statistic), serving as an asymptotically equivalent alternative
#' to the chi-squared test.
#'
#' \code{gTest} performs log-likelihood ratio contingency table tests and
#' goodness-of-fit tests.
#'
#' The G-test is also called "Likelihood Ratio Test" and is asymptotically
#' equivalent to the Pearson chi-squared test but not usually used when
#' analyzing 2x2 tables. It is used in logistic regression and loglinear
#' modeling which involves contingency tables.
#'
#' If \code{x} is a matrix with one row or column, or if \code{x} is a
#' vector and \code{y} is not given, then a \emph{goodness-of-fit test} is
#' performed (\code{x} is treated as a one-dimensional contingency table).
#' The entries of \code{x} must be non-negative integers. In this case, the
#' hypothesis tested is whether the population probabilities equal those in
#' \code{p}, or are all equal if \code{p} is not given.
#'
#' If \code{x} is a matrix with at least two rows and columns, it is taken
#' as a two-dimensional contingency table: the entries of \code{x} must be
#' non-negative integers. Otherwise, \code{x} and \code{y} must be vectors
#' or factors of the same length; cases with missing values are removed,
#' the objects are coerced to factors, and the contingency table is
#' computed from these. Then the G-test is performed on the null hypothesis
#' that the joint distribution of the cell counts in a 2-dimensional
#' contingency table is the product of the row and column marginals.
#'
#' Williams' correction (Williams, 1976) divides the statistic by a factor
#' \eqn{q > 1} and can be used for both test types. Yates' continuity
#' correction is only defined for 2x2 tables (independence) resp. two data
#' values (goodness-of-fit).
#'
#' @param x a numeric vector or matrix. \code{x} and \code{y} can also both
#' be factors.
#' @param y a numeric vector; ignored if \code{x} is a matrix. If \code{x}
#' is a factor, \code{y} should be a factor of the same length.
#' @param correct the correction to be applied, one of \code{"none"}
#' (default), \code{"williams"} or \code{"yates"}. See the Details.
#' @param p a vector of probabilities of the same length as \code{x}
#' (goodness-of-fit test only). An error is given if any entry of \code{p}
#' is negative.
#' @param rescaleP logical; if \code{TRUE} then \code{p} is rescaled (if
#' necessary) to sum to 1. If \code{rescaleP} is \code{FALSE}, and \code{p}
#' does not sum to 1, an error is given.
#' @return A list with class \code{"htest"} containing the following
#' components:
#' \item{statistic}{the value of the G test statistic.}
#' \item{parameter}{the degrees of freedom of the approximate chi-squared
#' distribution of the test statistic.}
#' \item{p.value}{the p-value of the test.}
#' \item{method}{a character string indicating the type of test performed,
#' and whether a correction was used.}
#' \item{data.name}{a character string giving the name(s) of the data.}
#' \item{observed}{the observed counts (before any continuity correction).}
#' \item{expected}{the expected counts under the null hypothesis.}
#'
#' @note Based on code by Pete Hurd, adapted to conform to package standards.
#'
#' @seealso \code{\link{chisq.test}}
#'
#' @references Agresti, A. (2007) \emph{An Introduction to Categorical Data
#' Analysis}, 2nd ed., New York: John Wiley & Sons. Page 38.
#'
#' Sokal, R. R. and Rohlf, F. J. (2012) \emph{Biometry: The Principles and
#' Practice of Statistics in Biological Research}, 4th ed., New York:
#' W. H. Freeman and Co.
#'
#' Williams, D. A. (1976) Improved likelihood ratio tests for complete
#' contingency tables. \emph{Biometrika}, 63, 33-37.
#'
#' @examples
#' ## From Agresti (2007), p. 39
#' M <- as.table(rbind(c(762, 327, 468), c(484, 239, 477)))
#' dimnames(M) <- list(gender = c("M", "F"),
#'                     party  = c("Democrat", "Independent", "Republican"))
#'
#' (Xsq <- gTest(M))   # Prints test summary
#'
#' Xsq$observed        # observed counts (same as M)
#' Xsq$expected        # expected counts under the null
#'
#'
#' ## Testing for population probabilities
#' ## Case A. Tabulated data
#' x <- c(A = 20, B = 15, C = 25)
#' gTest(x)
#' gTest(as.table(x))             # the same
#' x <- c(89, 37, 30, 28, 2)
#' p <- c(40, 20, 20, 15, 5)
#' try(
#' gTest(x, p = p)                # gives an error
#' )
#' # works
#' p <- c(0.40, 0.20, 0.20, 0.19, 0.01)
#' # Expected count in category 5
#' # is 1.86 < 5 ==> chi square approx.
#' gTest(x, p = p)                # maybe doubtful, but is ok!
#'
#' ## Case B. Raw data
#' x <- trunc(5 * runif(100))
#' gTest(table(x))                # NOT 'gTest(x)'!
#'
#' @family test.categorical
#' @concept categorical-test
#' @concept goodness-of-fit
#'
#' @export
gTest <- function(x, y = NULL, correct = c("none", "williams", "yates"),
                  p = rep(1 / length(x), length(x)), rescaleP = FALSE) {

  correct <- match.arg(correct)

  DNAME <- deparse1(substitute(x))

  if (is.data.frame(x)) x <- as.matrix(x)
  if (is.matrix(x)) {
    if (min(dim(x)) == 1L)
      x <- as.vector(x)
  }

  if (!is.matrix(x) && !is.null(y)) {
    if (length(x) != length(y))
      stop("'x' and 'y' must have the same length")
    DNAME <- paste(DNAME, "and", deparse1(substitute(y)))
    OK <- complete.cases(x, y)
    x <- as.factor(x[OK])
    y <- as.factor(y[OK])
    if (nlevels(x) < 2L || nlevels(y) < 2L)
      stop("'x' and 'y' must have at least 2 levels")
    x <- table(x, y)
  }

  if (any(x < 0) || anyNA(x))
    stop("all entries of 'x' must be nonnegative and finite")
  if ((n <- sum(x)) == 0)
    stop("at least one entry of 'x' must be positive")

  # the observed counts, before any continuity correction
  OBSERVED <- x

  if (is.matrix(x)) {

    ## ---- test of independence -----------------------------------------

    if (correct == "yates") {
      if (nrow(x) != 2L || ncol(x) != 2L)
        stop("Yates' correction requires a 2 x 2 matrix")
      # shift all cells towards their expected values by 0.5
      # (margins are preserved)
      if (x[1, 1] * x[2, 2] - x[1, 2] * x[2, 1] > 0) {
        x <- x + 0.5
        diag(x) <- diag(x) - 1
      } else {
        x <- x - 0.5
        diag(x) <- diag(x) + 1
      }
    }

    sr <- rowSums(x)
    sc <- colSums(x)
    E  <- outer(sr, sc) / n

    g <- sum(x[x > 0] * log(x[x > 0] / E[x > 0]))

    q <- if (correct == "williams")
      1 + ((n * sum(1 / sr) - 1) * (n * sum(1 / sc) - 1)) /
        (6 * n * (ncol(x) - 1) * (nrow(x) - 1))
    else 1

    STATISTIC <- 2 * g / q
    PARAMETER <- (nrow(x) - 1L) * (ncol(x) - 1L)

    METHOD <- paste("Log likelihood ratio (G-test) test of independence",
                    switch(correct,
                           none     = "without correction",
                           williams = "with Williams' correction",
                           yates    = "with Yates' correction"))

  } else {

    ## ---- goodness-of-fit test -------------------------------------------

    METHOD <- "Log likelihood ratio (G-test) goodness of fit test"

    if (length(dim(x)) > 2L)
      stop("invalid 'x'")
    if (length(x) == 1L)
      stop("'x' must at least have 2 elements")
    if (length(x) != length(p))
      stop("'x' and 'p' must have the same number of elements")
    if (any(p < 0))
      stop("probabilities must be non-negative")
    if (abs(sum(p) - 1) > sqrt(.Machine$double.eps)) {
      if (rescaleP)
        p <- p / sum(p)
      else
        stop("probabilities must sum to 1")
    }

    E <- n * p
    names(E) <- names(x)

    if (correct == "yates") {
      if (length(x) != 2L)
        stop("Yates' correction requires 2 data values")
      if (x[1] - E[1] > 0.25) {
        x[1] <- x[1] - 0.5
        x[2] <- x[2] + 0.5
      } else if (E[1] - x[1] > 0.25) {
        x[1] <- x[1] + 0.5
        x[2] <- x[2] - 0.5
      }
    }

    g <- sum(x[x > 0] * log(x[x > 0] / E[x > 0]))

    q <- if (correct == "williams")
      1 + (length(x) + 1) / (6 * n)
    else 1

    STATISTIC <- 2 * g / q
    PARAMETER <- length(x) - 1L
  }

  structure(
    list(statistic = c(G = STATISTIC),
         parameter = c(df = PARAMETER),
         p.value   = pchisq(STATISTIC, PARAMETER, lower.tail = FALSE),
         method    = METHOD,
         data.name = DNAME,
         observed  = OBSERVED,
         expected  = E),
    class = "htest")
}
