
#' Random proportions summing to 1
#'
#' Generates random positive values that sum to 1 by sampling from a
#' Dirichlet distribution with all concentration parameters equal to 1,
#' i.e. Dirichlet(1, ..., 1). This corresponds to a uniform distribution
#' over the simplex.
#'
#' Optionally, the result can be rounded to a specified number of digits
#' while preserving the total sum of 1.
#'
#' @param n Integer. Number of components (must be >= 1).
#' @param digits Optional integer. If provided, the result is rounded to
#'   the specified number of decimal places. A correction is applied to
#'   ensure the sum remains exactly 1.
#'
#' @return Numeric vector of length `n` with non-negative entries summing to 1.
#'
#' @details
#' Values are generated via normalized Gamma draws:
#' \deqn{x_i = g_i / \sum_j g_j, \quad g_i \sim \Gamma(1, 1)}
#'
#' If `digits` is specified, rounding may introduce small numerical
#' deviations. These are corrected by adjusting the largest component,
#' which minimizes relative distortion and avoids negative values.
#'
#' @examples
#' rsum1(5)
#' rsum1(5, digits = 2)
#' rsum1(1)
#'


#' @export
rsum1 <- function(n, digits = NULL) {
  x <- rdirichlet(1, rep(1, n))
  x <- as.numeric(x)
  
  if (!is.null(digits)) {
    x <- round(x, digits)
    x[which.max(x)] <- x[which.max(x)] + (1 - sum(x))
  }
  
  x
}

