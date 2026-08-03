
#' Cochran-Armitage Test for Detecting Trends in Proportions Across Ordered Groups
#'
#' A test for linear trend in proportions across ordered categories in
#' \eqn{2 \times k}{2xk} contingency tables, typically used in
#' epidemiological dose-response analyses.
#'
#' Perform a Cochran-Armitage test for trend in binomial proportions across
#' the levels of a single variable. This test is appropriate only when one
#' variable has two levels and the other variable is ordinal. The two-level
#' variable represents the response, and the other represents an
#' explanatory variable with ordered levels. The null hypothesis is the
#' hypothesis of no trend, which means that the binomial proportion is the
#' same for all levels of the explanatory variable.
#'
#' The table is oriented such that the ordered explanatory variable is in
#' the rows and the binary response in the columns (a \eqn{2 \times k}{2xk}
#' table is transposed accordingly). Row scores are taken from the row
#' dimnames if these are numeric, and are sequential integers otherwise.
#'
#' The alternatives refer to the trend in the proportion of the
#' \emph{first} response column: \code{"increasing"} tests whether this
#' proportion increases with the ordered levels, \code{"decreasing"} tests
#' the reverse.
#'
#' @param x a \eqn{k \times 2}{kx2} or \eqn{2 \times k}{2xk} frequency
#' table or matrix with nonnegative counts.
#' @param alternative a character string specifying the alternative
#' hypothesis, must be one of \code{"two.sided"} (default),
#' \code{"increasing"} or \code{"decreasing"}. You can specify just the
#' initial letter. See the Details for the direction convention.
#'
#' @return A list of class \code{"htest"}, containing the following
#' components:
#'   \item{\code{statistic}}{the z-statistic of the test.}
#'   \item{\code{parameter}}{the number of levels of the ordered variable
#'     (named \code{dim}).}
#'   \item{\code{p.value}}{the p-value for the test.}
#'   \item{\code{alternative}}{a character string describing the alternative
#'     hypothesis.}
#'   \item{\code{method}}{the character string "Cochran-Armitage test for
#'     trend".}
#'   \item{\code{data.name}}{a character string giving the names of the data.}
#'
#' @note
#' Based on code by Eric Lecoutre, adapted to conform to package standards.
#'
#' \url{https://stat.ethz.ch/pipermail/r-help/2005-July/076371.html}
#'
#' Results are consistent with SAS PROC FREQ. They may differ slightly from
#' \pkg{coin}'s \code{independence_test(..., teststat = "scalar")}, which
#' uses a different variance formula.
#'
#' @seealso
#' \code{\link{prop.trend.test}},
#' \href{https://documentation.sas.com/doc/en/pgmsascdc/9.4_3.5/procstat/procstat_freq_details76.htm}{SAS PROC FREQ documentation}
#'
#' @references Agresti, A. (2002) \emph{Categorical Data Analysis}. John
#' Wiley & Sons.
#'
#' @examples
#' # http://www.lexjansen.com/pharmasug/2007/sp/sp05.pdf, pp. 4
#' dose <- matrix(c(10,9,10,7, 0,1,0,3), byrow=TRUE, nrow=2,
#'                dimnames=list(resp=0:1, dose=0:3))
#'
#' cochranArmitageTest(dose)
#' cochranArmitageTest(dose, alternative = "increasing")
#'
#' # compare with coin::independence_test(..., teststat = "scalar")
#' # (see the Note on the variance formula)
#' lungtumor <- data.frame(dose = rep(c(0, 1, 2), c(40, 50, 48)),
#'                         tumor = c(rep(c(0, 1), c(38, 2)),
#'                                   rep(c(0, 1), c(43, 7)),
#'                                   rep(c(0, 1), c(33, 15))))
#' tab <- table(lungtumor$dose, lungtumor$tumor)
#' cochranArmitageTest(tab)
#'
#' # similar to prop.trend.test (uses integer scores 1..k instead of dimnames)
#' prop.trend.test(tab[,1], apply(tab, 1, sum))
#'
#' # SAS PROC FREQ reference
#' # https://support.sas.com/documentation/onlinedoc/stat/142/freq.pdf, pp 2868
#' pain <- structure(c(26, 6, 26, 7, 23, 9, 18, 14, 9, 23),
#'                   dim = c(2L, 5L),
#'                   dimnames = list(adverse = c("No", "Yes"),
#'                                   dose    = c("0", "1", "2", "3", "4")),
#'                   class = "table")
#'
#' cochranArmitageTest(pain)
#'
#' @family test.trend
#' @concept categorical-test
#' @concept nonparametric
#'
#' @export
cochranArmitageTest <- function(x,
                                alternative = c("two.sided", "increasing",
                                                "decreasing")) {

  alternative <- match.arg(alternative)

  DNAME <- deparse1(substitute(x))

  x <- as.matrix(x)

  if (length(dim(x)) != 2L || !any(dim(x) == 2L))
    stop("Cochran-Armitage test for trend must be used with a kx2 table",
         call. = FALSE)
  if (anyNA(x) || any(x < 0))
    stop("all entries of 'x' must be nonnegative counts", call. = FALSE)

  # ensure rows are the ordered variable (k levels),
  # columns are the response (2 levels)
  if (dim(x)[2L] != 2L)
    x <- t(x)

  nidot <- apply(x, 1L, sum)
  n     <- sum(nidot)

  # row scores from dimnames if numeric, otherwise sequential integers
  Ri   <- scores(x, MARGIN = 1L, "table")
  Rbar <- sum(nidot * Ri) / n
  s2   <- sum(nidot * (Ri - Rbar)^2)

  pdot1 <- sum(x[, 1L]) / n
  z     <- sum(x[, 1L] * (Ri - Rbar)) / sqrt(pdot1 * (1 - pdot1) * s2)

  PVAL <- switch(
    alternative,
    "two.sided"  = 2 * pnorm(abs(z), lower.tail = FALSE),
    "increasing" = pnorm(z, lower.tail = FALSE),
    "decreasing" = pnorm(z)
  )

  structure(
    list(
      statistic   = c(Z = z),
      parameter   = c(dim = dim(x)[1L]),
      p.value     = PVAL,
      alternative = alternative,
      method      = "Cochran-Armitage test for trend",
      data.name   = DNAME
    ),
    class = "htest"
  )
}
