
#' Cramer-Von Mises Test for Assessing Deviations From Normality
#'
#' A goodness-of-fit test based on the integrated squared discrepancies
#' between the empirical and theoretical distribution functions. Similar
#' to the Anderson-Darling test, but without increased weighting of the
#' distribution tails.
#'
#' Performs the Cramer-von Mises test for the composite hypothesis of
#' normality, see e.g. Thode (2002, Sec. 5.1.3).
#'
#' The Cramer-von Mises test is an EDF omnibus test for the composite
#' hypothesis of normality. The test statistic is
#' \deqn{W = \frac{1}{12 n} + \sum_{i=1}^{n} \left(p_{(i)} - \frac{2i - 1}{2n}\right)^2,}{W = 1/(12n) + sum(p_(i) - (2i-1)/(2n))^2,}
#' where \eqn{p_{(i)} = \Phi([x_{(i)} - \bar{x}]/s)}. Here, \eqn{\Phi} is
#' the cumulative distribution function of the standard normal
#' distribution, and \eqn{\bar{x}} and \eqn{s} are mean and standard
#' deviation of the data values. The p-value is computed from the modified
#' statistic \eqn{Z = W (1 + 0.5/n)} according to Table 4.9 in Stephens
#' (1986).
#'
#' Missing values are silently removed.
#'
#' @param x a numeric vector of data values, the number of which must be
#' at least 8.
#' @return A list with class \code{"htest"} containing the following
#' components:
#'   \item{\code{statistic}}{the value of the Cramer-von Mises statistic.}
#'   \item{\code{p.value}}{the p-value of the test.}
#'   \item{\code{method}}{a character string indicating the test performed.}
#'   \item{\code{data.name}}{a character string giving the name of the data.}
#'
#' @note
#' Based on code by Juergen Gross previously published in the
#' \pkg{nortest} package, adapted to conform to package standards.
#'
#' @references Stephens, M.A. (1986) Tests based on EDF statistics. In:
#' D'Agostino, R.B. and Stephens, M.A., eds.: \emph{Goodness-of-Fit
#' Techniques}. New York: Marcel Dekker.
#'
#' Thode Jr., H.C. (2002) \emph{Testing for Normality}. New York: Marcel
#' Dekker.
#'
#' @seealso \code{\link{shapiro.test}} for performing the Shapiro-Wilk test
#' for normality, \code{\link{andersonDarlingTest}},
#' \code{pharos::plotQQ()} for producing extended normal quantile-quantile
#' plots
#'
#' @examples
#' set.seed(1)
#' cramerVonMisesTest(rnorm(100, mean = 5, sd = 3))
#' cramerVonMisesTest(runif(100, min = 2, max = 4))
#'
#' @family test.normality
#' @concept normality-test
#' @concept goodness-of-fit
#'
#' @export
cramerVonMisesTest <- function(x) {

  DNAME <- deparse1(substitute(x))

  x <- sort(x[complete.cases(x)])
  n <- length(x)

  if (n < 8)
    stop("sample size must be at least 8")

  p <- pnorm((x - mean(x)) / sd(x))

  W <- 1 / (12 * n) +
    sum((p - (2 * seq_len(n) - 1) / (2 * n))^2)

  # modified statistic and p-value approximation,
  # Stephens (1986), Table 4.9
  WW <- (1 + 0.5 / n) * W

  if (WW < 0.0275) {
    pval <- 1 - exp(-13.953 + 775.5 * WW - 12542.61 * WW^2)
  } else if (WW < 0.051) {
    pval <- 1 - exp(-5.903 + 179.546 * WW - 1515.29 * WW^2)
  } else if (WW < 0.092) {
    pval <- exp(0.886 - 31.62 * WW + 10.897 * WW^2)
  } else if (WW < 1.1) {
    pval <- exp(1.111 - 34.242 * WW + 12.832 * WW^2)
  } else {
    warning("p-value is smaller than 7.37e-10, ",
            "cannot be computed more accurately")
    pval <- 7.37e-10
  }

  structure(
    list(statistic = c(W = W),
         p.value   = pval,
         method    = "Cramer-von Mises normality test",
         data.name = DNAME),
    class = "htest")
}
