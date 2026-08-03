
#' Lehmacher's Test for Locating Asymmetry in Paired Multicategory Data
#'
#' A nonparametric test for marginal homogeneity in square contingency
#' tables for dependent samples, based on a normal approximation of the
#' cell frequency differences.
#'
#' Performs Lehmacher's chi-squared test for marginal homogeneity in a
#' square two-dimensional contingency table.
#'
#' Unlike Bowker's test of symmetry, which tests whether
#' \eqn{P(i,j) = P(j,i)} for every off-diagonal cell pair, Lehmacher's
#' test addresses marginal homogeneity: the null hypothesis is that, for
#' every category \eqn{i}, the row and column marginal probabilities agree,
#' \eqn{P(i \cdot) = P(\cdot i)}. One test statistic is computed per
#' category and referred to a chi-squared distribution with 1 degree of
#' freedom; the resulting p-values are adjusted for multiple comparisons.
#'
#' If \code{x} is a matrix, it is taken as a two-dimensional contingency
#' table, and hence its entries should be nonnegative integers. Otherwise,
#' both \code{x} and \code{y} must be vectors or factors of the same
#' length. Incomplete cases are removed, vectors are coerced into factors,
#' and the contingency table is computed from these.
#'
#' @name lehmacherTest
#' @param x either a two-dimensional square contingency table in matrix
#' form, or a factor object.
#' @param y a factor object; ignored if \code{x} is a matrix.
#' @param p.adjust.method the method used to adjust the per-category
#' p-values for multiple comparisons, passed to \code{\link{p.adjust}}.
#' Default is \code{"hochberg"}, as recommended by Lehmacher (1980).
#' @return A list with class \code{c("MHTest", "htest")} containing the
#' following components:
#' \item{statistic}{a vector with the value of the test statistic for each
#' category.}
#' \item{parameter}{the degrees of freedom, which is always 1.}
#' \item{p.value}{a vector with the p-values of the individual tests.}
#' \item{p.value.corr}{a vector with the adjusted p-values of the
#' individual tests (see \code{p.adjust.method}).}
#' \item{method}{a character string indicating the test performed.}
#' \item{data.name}{a character string giving the name of the data.}
#'
#' @seealso \code{\link{stuartMaxwellTest}}, \code{\link{bhapkarTest}} for
#' related tests of marginal homogeneity; \code{\link{mcnemar.test}} for
#' the 2x2 case
#'
#' @references Lehmacher, W. (1980) Simultaneous sign tests for marginal
#' homogeneity of square contingency tables. \emph{Biometrical Journal},
#' 22 (8), 795-798.
#'
#' @examples
#' x <- matrix(c(400, 40, 20, 10,
#'               50, 300, 60, 20,
#'               10, 40, 120, 5,
#'               5, 90, 50, 80), nrow = 4, byrow = TRUE,
#'             dimnames = list(LETTERS[1:4], LETTERS[1:4]))
#'
#' lehmacherTest(x)
#'
#' @family test.categorical
#' @concept categorical-test
#' @concept marginal-homogeneity
#' @concept paired-samples
#'
#' @export
lehmacherTest <- function(x, y = NULL, p.adjust.method = "hochberg") {

  CT <- resolveContingency(
    x = x,
    y = y,
    square = TRUE
  )

  x <- CT$table
  DNAME <- CT$data.name

  rsum <- rowSums(x)
  csum <- colSums(x)

  denom <- rsum + csum - 2 * diag(x)

  STATISTIC <- ifelse(
    denom > 0,
    (rsum - csum)^2 / denom,
    0
  )

  PVAL <- pchisq(
    STATISTIC,
    df = 1L,
    lower.tail = FALSE
  )

  structure(
    list(
      statistic    = STATISTIC,
      parameter    = c(df = 1L),
      p.value      = PVAL,
      p.value.corr = p.adjust(PVAL, method = p.adjust.method),
      method       = "Lehmacher test for marginal homogeneity",
      data.name    = DNAME
    ),
    class = c("MHTest", "htest")
  )
}



#' @param digits a non-null value for digits specifies the minimum number
#' of significant digits to be printed. See \code{\link{print.default}}.
#' @param \dots further arguments to be passed to or from other methods,
#' ignored in this function.
#'
#' @rdname lehmacherTest
#' @export
print.MHTest <- function(x, digits = 1L, ...) {

  cat("\n")
  cat(strwrap(x$method, prefix = "\t"), sep = "\n")
  cat("\n")
  cat("data:  ", x$data.name, "\n\n", sep = "")

  out <- cbind(
    fm(x$statistic, digits = digits),
    fm(x$p.value, fmt = "p"),
    fm(x$p.value.corr, fmt = "p"),
    fm(x$p.value.corr, fmt = "*")
  )

  colnames(out) <- c("Chi\u00B2", "p-value", "p-adj", " ")
  rownames(out) <- if (!is.null(names(x$statistic)))
    names(x$statistic)
  else
    seq_along(x$statistic)

  print.default(out, digits = 3L, print.gap = 3, quote = FALSE, right = TRUE)
  .printSignifCodes()
  cat("\n")

  invisible(x)
}
