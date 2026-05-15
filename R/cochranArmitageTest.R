
#' Cochran-Armitage Test for Trend
#'
#' A test for linear trend in proportions across ordered categories in 2×k
#' contingency tables, typically used in epidemiological dose-response analyses.
#'
#' Perform a Cochran Armitage test for trend in binomial proportions across the
#' levels of a single variable. This test is appropriate only when one variable
#' has two levels and the other variable is ordinal. The two-level variable
#' represents the response, and the other represents an explanatory variable
#' with ordered levels. The null hypothesis is the hypothesis of no trend,
#' which means that the binomial proportion is the same for all levels of the
#' explanatory variable.
#'
#'
#' @param x a frequency table or a matrix.
#' @param alternative a character string specifying the alternative hypothesis,
#' must be one of \code{"two.sided"} (default), \code{"increasing"} or
#' \code{"decreasing"}. You can specify just the initial letter.
#'
#' @return A list of class \code{htest}, containing the following components:
#' \item{statistic}{ the z-statistic of the test.}
#' \item{parameter}{ the dimension of the table.}
#' \item{p.value}{ the p-value for the test.}
#' \item{alternative}{a character string describing the alternative hypothesis.}
#' \item{method}{the character string \dQuote{Cochran-Armitage test for trend}.}
#' \item{data.name}{a character string giving the names of the data.}
#'
#' @note
#' Based on code by Eric Lecoutre.
#'
#' \url{https://stat.ethz.ch/pipermail/r-help/2005-July/076371.html}
#'
#' Results are consistent with SAS PROC FREQ. They may differ from
#' \pkg{coin}'s \code{independence_test(..., teststat = "scalar")},
#' which uses a different variance formula.
#' 
#' @seealso
#' \code{\link{prop.trend.test}}
#' \href{https://documentation.sas.com/doc/en/pgmsascdc/9.4_3.5/procstat/procstat_freq_details76.htm}{SAS PROC FREQ documentation}
#'
#' @references Agresti, A. (2002) \emph{Categorical Data Analysis}. John Wiley
#' & Sons
#'
#' @examples
#'
#' # http://www.lexjansen.com/pharmasug/2007/sp/sp05.pdf, pp. 4
#' dose <- matrix(c(10,9,10,7, 0,1,0,3), byrow=TRUE, nrow=2,
#'                dimnames=list(resp=0:1, dose=0:3))
#'
#' cochranArmitageTest(dose)
#' cochranArmitageTest(dose, alternative = "increasing")
#'
#' # consistent with coin::independence_test(..., teststat = "scalar")
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
#' # SAS PROC FREQ reference (note: SAS uses a different variance formula)
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
#' @concept hypothesis-testing
#' @concept nonparametric
#' @concept table-manipulation
#'


#' @export
cochranArmitageTest <- function(x,
                                alternative = c("two.sided", "increasing", "decreasing")) {
  
  DNAME <- deparse(substitute(x))
  
  if (!any(dim(x) == 2L))
    stop("Cochran-Armitage test for trend must be used with rx2-table",
         call. = FALSE)
  
  # Ensure rows are the ordered variable (k levels), columns are the response (2 levels)
  if (dim(x)[2L] != 2L)
    x <- t(x)
  
  nidot <- apply(x, 1L, sum)
  n     <- sum(nidot)
  
  # Row scores from dimnames if numeric, otherwise sequential integers
  Ri   <- scores(x, MARGIN = 1L, "table")
  Rbar <- sum(nidot * Ri) / n
  s2   <- sum(nidot * (Ri - Rbar)^2)
  
  pdot1 <- sum(x[, 1L]) / n
  z     <- sum(x[, 1L] * (Ri - Rbar)) / sqrt(pdot1 * (1 - pdot1) * s2)
  
  alternative <- match.arg(alternative)
  
  PVAL <- switch(
    alternative,
    "two.sided"  = 2 * pnorm(abs(z), lower.tail = FALSE),
    "increasing" = pnorm(z, lower.tail = FALSE),
    "decreasing" = pnorm(z)
  )
  
  STATISTIC <- c(Z = z)
  PARAMETER <- c(dim = dim(x)[1L])
  
  structure(
    list(
      statistic   = STATISTIC,
      parameter   = PARAMETER,
      p.value     = PVAL,
      alternative = alternative,
      method      = "Cochran-Armitage test for trend",
      data.name   = DNAME
    ),
    class = "htest"
  )
}

