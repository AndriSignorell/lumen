
#' Fast correlation test for a matrix
#'
#' A test for significant correlation between two variables, supporting 
#' Pearson, Spearman, and Kendall correlation coefficients.
#' 
#' Computes Pearson correlations, pairwise sample sizes, and p-values for all
#' variable pairs in a numeric matrix. The p-values are computed directly from
#' the correlation coefficients using the exact Beta distribution under the
#' null hypothesis of zero correlation.
#'
#' Compared to repeatedly calling `stats::cor.test()`, this implementation is
#' fully vectorised and substantially faster for matrices with many variables.
#'
#' Missing values are handled using pairwise complete observations.
#'
#' @param X Numeric matrix or data.frame.
#' @param use Character string passed to `stats::cor()`. Default is
#'   `"pairwise.complete.obs"`.
#' @param triangle Which part of the matrix should be returned:
#'   `"full"` (default), `"upper"`, or `"lower"`.
#' @param sig.level Optional significance threshold. If supplied,
#'   correlations with p-values larger than `sig.level` are replaced with `NA`.
#'
#' @return A list with three matrices:
#' \describe{
#'   \item{cor}{Correlation matrix}
#'   \item{p}{Matrix of p-values}
#'   \item{n}{Matrix of pairwise sample sizes}
#' }
#'
#' @details
#' For a Pearson correlation coefficient \(r\) with sample size \(n\),
#' the two-sided p-value under the null hypothesis rho = 0 can be computed as
#'
#' \deqn{p = 2 * F[Beta]((1 - |r|)/2, 1/2, (n-2)/2)}
#'
#' where F[Beta] is the Beta distribution CDF.
#' 
#' This formulation avoids explicitly computing the t-statistic and is
#' numerically stable for correlations close to ±1.
#'
#' @examples
#' X <- matrix(rnorm(200), 50, 4)
#' res <- corTest(X)
#'
#' # Only significant correlations in the upper triangle
#' res2 <- corTest(X, triangle = "upper", sig.level = 0.05)
#'
#' @family test.correlation
#' @concept hypothesis-testing
#' @concept correlation
#' @concept descriptive-statistics
#'
#'
#' @export
corTest <- function(X,
                          use = "pairwise.complete.obs",
                          triangle = c("full", "upper", "lower"),
                          sig.level = NULL) {
  

  # approach for tri matrix only in oneliner:
  #
  # naIf(m * (row(m) <= col(m)), 0)
  # {m[lower.tri(m)] <- NA; m}
  
  
  triangle <- match.arg(triangle)
  
  X <- as.matrix(X)
  
  R <- cor(X, use = use)
  
  n <- crossprod(!is.na(X))
  
  P <- 2 * pbeta((1 - abs(R)) / 2, 0.5, (n - 2) / 2)
  
  diag(P) <- NA
  
  if (!is.null(sig.level)) {
    R[P > sig.level] <- NA
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
    cor = R,
    p = P,
    n = n
  )
}
