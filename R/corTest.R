
#' Fast Correlation Test for a Matrix
#'
#' Compute correlation coefficients, p-values, and pairwise sample sizes
#' for all variable pairs in a numeric matrix.
#'
#' Pearson, Spearman, and Kendall correlations are supported via the
#' \code{method} argument passed to \code{\link{cor}}.
#'
#' Compared to repeatedly calling \code{\link[stats]{cor.test}},
#' this implementation is fully vectorised and substantially faster
#' for matrices with many variables.
#'
#' Missing values are handled using pairwise complete observations.
#'
#' @param x Numeric matrix or data.frame.
#' @param method Character string specifying the correlation method.
#'   Passed to \code{\link{cor}}.
#'
#'   One of:
#'
#'   \describe{
#'     \item{\code{\"pearson\"}}{Pearson product-moment correlation}
#'     \item{\code{\"spearman\"}}{Spearman rank correlation}
#'     \item{\code{\"kendall\"}}{Kendall rank correlation}
#'   }
#'
#'   Default is \code{\"pearson\"}.
#'
#' @param use Character string passed to \code{\link{cor}}.
#'   Default is \code{\"pairwise.complete.obs\"}.
#'
#' @param triangle Character string specifying which part of the
#'   matrices should be returned.
#'
#'   One of:
#'
#'   \describe{
#'     \item{\code{\"full\"}}{Return the complete matrices}
#'     \item{\code{\"upper\"}}{Return only the upper triangle}
#'     \item{\code{\"lower\"}}{Return only the lower triangle}
#'   }
#'
#'   Default is \code{\"full\"}.
#'
#' @param maxPValue Optional upper limit for displayed p-values.
#'   Correlations with p-values larger than \code{maxPValue}
#'   are replaced with \code{NA}.
#'
#' @return
#' A list with three matrices:
#'
#' \describe{
#'   \item{cor}{Correlation matrix}
#'   \item{pValue}{Matrix of p-values}
#'   \item{n}{Matrix of pairwise sample sizes}
#' }
#'
#' @details
#' For a Pearson correlation coefficient \eqn{r}
#' with sample size \eqn{n},
#' the two-sided p-value under the null hypothesis
#' \eqn{\rho = 0} can be computed as:
#'
#' \deqn{
#' p = 2 * F_{Beta}((1 - |r|)/2, 1/2, (n-2)/2)
#' }
#'
#' where \eqn{F_{Beta}} is the Beta distribution CDF.
#'
#' This formulation avoids explicitly computing the t-statistic and is
#' numerically stable for correlations close to ±1.
#'
#' @examples
#' X <- matrix(rnorm(200), 50, 4)
#'
#' res <- corTest(X)
#'
#' # only significant correlations in the upper triangle
#' res2 <- corTest(
#'   X,
#'   triangle = "upper",
#'   maxPValue = 0.05
#' )
#'

#' @family test.correlation  
#' @concept correlation  
#' @concept hypothesis-test
#'
#'
#' @export
corTest <- function(
    x,
    method = c("pearson", "spearman", "kendall"),
    use = "pairwise.complete.obs",
    triangle = c("full", "upper", "lower"),
    maxPValue = NULL
) {
  
  # approach for triangle matrix only in one line:
  #
  # naIf(m * (row(m) <= col(m)), 0)
  # {m[lower.tri(m)] <- NA; m}
  
  method   <- match.arg(method)
  triangle <- match.arg(triangle)
  
  x <- as.matrix(x)
  
  R <- cor(x, method = method, use = use)
  
  n <- crossprod(!is.na(x))
  
  # Correct two-sided p-value via t-statistic: t = r*sqrt((n-2)/(1-r^2))
  # p = 2 * pt(-|t|, df = n-2)
  # Using pt(x, df) = pbeta(df/(df+x^2), df/2, 1/2)/2 for the lower tail:
  # 2 * pt(-|t|, n-2) = pbeta((n-2) / ((n-2) + (n-2)*R^2/(1-R^2)), (n-2)/2, 1/2)
  # Simpler: compute T matrix directly
  T_stat <- R * sqrt((n - 2) / pmax(1 - R^2, .Machine$double.eps))
  P      <- 2 * pt(-abs(T_stat), df = n - 2)
  
  diag(P) <- NA
  
  if (!is.null(maxPValue)) {
    R[P > maxPValue] <- NA
  }
  
  if (triangle == "upper") {
    
    R[lower.tri(R)] <- NA
    P[lower.tri(P)] <- NA
    n[lower.tri(n)] <- NA
  }
  
  if (triangle == "lower") {
    
    R[upper.tri(R)] <- NA
    P[upper.tri(P)] <- NA
    n[upper.tri(n)] <- NA
  }
  
  list(
    cor    = R,
    pValue = P,
    n      = n
  )
}

