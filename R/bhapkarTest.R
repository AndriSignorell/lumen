
#' Bhapkar Marginal Homogeneity Test for Comparing Paired Marginal Distributions Across Multiple Categories
#'
#' An asymptotic chi-squared test for marginal homogeneity in square
#' contingency tables for dependent samples, similar to the Stuart-Maxwell
#' test but based on a different test statistic and generally slightly more
#' powerful.
#'
#' Bhapkar's test (Bhapkar, 1966) is used to assess marginal homogeneity
#' in square contingency tables. It is based on the asymptotic normality
#' of marginal proportions and is closely related to the generalized
#' McNemar test, as implemented in \code{\link{stuartMaxwellTest}}.
#'
#' The two tests differ only in the estimation of the variance-covariance
#' matrix of the marginal proportions and are asymptotically equivalent
#' (Keefe, 1982), meaning that for large sample sizes they yield the same
#' chi-squared statistic. For finite samples, however, the Bhapkar test is
#' generally more powerful and is therefore preferred in practice.
#'
#' @param x either a 2-way \eqn{k \times k}{k x k} contingency table in
#' matrix form, or a factor.
#' @param y a factor with the same levels as \code{x}; ignored if \code{x}
#' is a matrix.
#' @return A list with class \code{"htest"} containing the following
#' components:
#'   \item{\code{statistic}}{the value of the chi-squared test statistic.}
#'   \item{\code{parameter}}{the degrees of freedom of the approximate
#'     chi-squared distribution of the test statistic.}
#'   \item{\code{p.value}}{the p-value of the test.}
#'   \item{\code{method}}{a character string indicating the test performed.}
#'   \item{\code{data.name}}{a character string giving the name of the data.}
#'
#' @references Bhapkar V.P. (1966) A note on the equivalence of two test
#' criteria for hypotheses in categorical data. \emph{Journal of the
#' American Statistical Association}, 61: 228-235.
#'
#' Ireland C.T., Ku H.H., and Kullback S. (1969) Symmetry and marginal
#' homogeneity of an r x r contingency table. \emph{Journal of the American
#' Statistical Association}, 64: 1323-1341.
#'
#' Keefe T.J. (1982) On the relationship between two tests for homogeneity
#' of the marginal distributions in a two-way classification.
#' \emph{Biometrika}, 69: 683-684.
#'
#' Sun X., Yang Z. (2008) Generalized McNemar's Test for Homogeneity of the
#' Marginal Distributions. \emph{SAS Global Forum 2008: Statistics and Data
#' Analysis}, Paper 382-208.
#'
#' @seealso [mcnemar.test],[chisq.test]
#'
#' @examples
#' # Source: https://john-uebersax.com/stat/mcnemar.htm#bhapkar
#' mc <- as.table(matrix(c(20,3,0,10,30,5,5,15,40), nrow=3))
#'
#' bhapkarTest(mc)
#'
#' @family test.categorical
#' @concept categorical-test
#' @concept marginal-homogeneity
#' @concept paired-samples
#'
#' @export
bhapkarTest <- function(x, y = NULL) {

  DNAME <- if (is.null(y)) deparse(substitute(x)) else
    paste(deparse(substitute(x)), "and", deparse(substitute(y)))

  res <- stuartMaxwellTest(x = x, y = y)

  Q_SM <- unname(res$statistic)
  Q_B  <- Q_SM / (1 - Q_SM / res$n)

  res$statistic <- c("chi-squared" = Q_B)
  res$p.value   <- pchisq(Q_B, df = res$parameter, lower.tail = FALSE)
  res$method    <- "Bhapkar test for marginal homogeneity"
  res$data.name <- DNAME

  res
}
