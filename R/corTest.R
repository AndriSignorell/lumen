
#' Fast Correlation Test  for Testing Pairwise Correlations in a Matrix
#'
#' Compute correlation coefficients, p-values, and pairwise sample sizes
#' for all variable pairs in a numeric matrix.
#'
#' Pearson, Spearman, and Kendall correlations are supported via the
#' `method` argument passed to [cor()].
#'
#' Compared to repeatedly calling [cor.test()], this
#' implementation is fully vectorised and substantially faster for matrices
#' with many variables.
#'
#' For Pearson and Spearman correlations, the two-sided p-value under the
#' null hypothesis \eqn{\rho = 0} is computed from the t statistic
#' \deqn{t = r \sqrt{\frac{n - 2}{1 - r^2}}}{t = r * sqrt((n-2)/(1-r^2))}
#' with \eqn{n - 2} degrees of freedom. This is exact for the Pearson
#' coefficient under normality and an approximation for the Spearman
#' coefficient (as used by `cor.test` for larger samples). For the
#' Kendall coefficient the normal approximation
#' \deqn{z = \frac{3 \tau \sqrt{n (n-1)}}{\sqrt{2 (2 n + 5)}}}{z = 3*tau*sqrt(n(n-1)) / sqrt(2(2n+5))}
#' is used; note that no correction for ties is applied here. Cells with
#' fewer than 3 (Pearson/Spearman) resp. 2 (Kendall) pairwise observations
#' get an `NA` p-value.
#'
#' If `use` requests pairwise handling of missing values (the
#' default), the sample sizes in `n` are the pairwise counts;
#' otherwise all pairs share the number of complete cases.
#'
#' @param x a numeric matrix or data frame.
#' @param method a character string specifying the correlation method, one
#' of `"pearson"` (default), `"spearman"` or `"kendall"`.
#' Passed to [cor()].
#' @param use a character string giving a method for computing correlations
#' in the presence of missing values, passed to [cor()]. Default
#' is `"pairwise.complete.obs"`.
#' @param triangle a character string specifying which part of the matrices
#' should be returned, one of `"full"` (default), `"upper"` or
#' `"lower"`. The other triangle is set to `NA`.
#' @param maxPValue optional upper limit for reported correlations.
#' Correlations with p-values larger than `maxPValue` are replaced
#' with `NA`.
#' 
#' @return A list with three matrices:
#'   \item{`cor`}{the correlation matrix.}
#'   \item{`pValue`}{the matrix of two-sided p-values (`NA` on the
#'     diagonal).}
#'   \item{`n`}{the matrix of sample sizes used per pair.}
#'
#' @seealso [cor()], [cor.test()]
#'
#' @examples
#' set.seed(1)
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
#'
#' @export
corTest <- function(x,
                    method = c("pearson", "spearman", "kendall"),
                    use = "pairwise.complete.obs",
                    triangle = c("full", "upper", "lower"),
                    maxPValue = NULL) {

  method   <- match.arg(method)
  triangle <- match.arg(triangle)

  x <- as.matrix(x)

  R <- cor(x, method = method, use = use)

  # sample sizes, consistent with the missing-value handling in cor()
  n <- if (grepl("^pairwise", use)) {
    crossprod(!is.na(x))
  } else {
    matrix(sum(complete.cases(x)), ncol(x), ncol(x),
           dimnames = dimnames(R))
  }

  if (method == "kendall") {
    # normal approximation for Kendall's tau (no ties correction)
    Z <- 3 * R * sqrt(n * (n - 1)) / sqrt(2 * (2 * n + 5))
    P <- 2 * pnorm(-abs(Z))
    P[n < 2] <- NA_real_
  } else {
    # exact for Pearson (under normality), approximate for Spearman:
    # t = r * sqrt((n-2)/(1-r^2)) with n-2 degrees of freedom
    T_stat <- R * sqrt((n - 2) / pmax(1 - R^2, .Machine$double.eps))
    P <- 2 * pt(-abs(T_stat), df = n - 2)
    P[n < 3] <- NA_real_
  }

  diag(P) <- NA

  if (!is.null(maxPValue))
    R[P > maxPValue] <- NA

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
