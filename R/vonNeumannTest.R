
#' Von Neumann's Successive Difference Test
#'
#' A test for randomness or autocorrelation in a sequence, based on the mean
#' square successive difference relative to the sample variance, closely
#' related to the Durbin-Watson test.
#'
#' The test is based on the von Neumann ratio statistic.
#'
#' The VN test statistic is in the unbiased case
#' \deqn{VN=\frac{\sum_{i=1}^{n-1}(x_i-x_{i+1})^2 \cdot
#' n}{\sum_{i=1}^{n}\left(x_i-\bar{x}\right)^2 \cdot (n-1)}}{
#' VN = sum(x_i - x_{i+1})^2 * n / (sum(x_i - mean(x))^2 * (n-1))}
#'
#' It is known that \eqn{(VN-\mu)/\sigma} is asymptotically standard normal,
#' where \eqn{\mu = 2n/(n-1)} and
#' \eqn{\sigma^2 = 4 n^2 (n-2) / [(n+1)(n-1)^3]}.
#'
#' The VN test statistic is in the original (biased) case
#' \deqn{VN=\frac{\sum_{i=1}^{n-1}(x_i-x_{i+1})^2}{
#' \sum_{i=1}^{n}\left(x_i-\bar{x}\right)^2}}{
#' VN = sum(x_i - x_{i+1})^2 / sum(x_i - mean(x))^2}
#'
#' The test statistic \eqn{(VN-2)/\sigma} is asymptotically standard normal,
#' where \eqn{\sigma^2 = 4(n-2) / [(n+1)(n-1)]}.
#'
#' Missing values are silently removed.
#'
#' @name vonNeumannTest
#'
#' @param x a numeric vector containing the observations.
#' @param alternative a character string specifying the alternative hypothesis,
#'   must be one of \code{"two.sided"} (default), \code{"greater"} or
#'   \code{"less"}.
#' @param unbiased logical. If \code{TRUE} (default), applies the
#'   finite-sample correction \eqn{n/(n-1)} so that VN is an unbiased
#'   estimate of the population value.
#'   
#' @return A list with class \code{"htest"} containing:
#' \item{statistic}{the normalized z-statistic.}
#' \item{parameter}{named vector with \code{n}, the sample size after
#'   removal of \code{NA}s.}
#' \item{p.value}{the p-value of the test.}
#' \item{alternative}{a character string describing the alternative
#'   hypothesis.}
#' \item{data.name}{a character string giving the name of the data.}
#' \item{vn}{the value of the VN statistic (not printed).}
#'
#' @seealso \code{\link{bartelsRankTest}}
#'
#' @references
#' von Neumann, J. (1941) Distribution of the ratio of the mean square
#' successive difference to the variance. \emph{Annals of Mathematical
#' Statistics} \bold{12}, 367--395.
#'
#' Young, L. C. (1941) On randomness in ordered sequences.
#' \emph{Annals of Mathematical Statistics} \bold{12}, 293--300.
#'
#' Bartels, R. (1982) The Rank Version of von Neumann's Ratio Test for
#' Randomness. \emph{Journal of the American Statistical Association},
#' \bold{77}(377), 40--46.
#'
#' @examples
#' set.seed(2)
#' vonNeumannTest(runif(20))
#'
#' # trend: small VN expected
#' vonNeumannTest(cumsum(rnorm(30)), alternative = "less")
#'
#' @rdname vonNeumannTest
#' @family test.correlation
#' @concept hypothesis-testing
#' @concept nonparametric
#' @concept time-series


#' @export
vonNeumannTest <- function(x,
                           alternative = c("two.sided", "less", "greater"),
                           unbiased    = TRUE) {
  
  DNAME <- deparse(substitute(x))
  
  x <- x[!is.na(x)]
  
  if (!is.numeric(x))
    stop("'x' must be a numeric vector")
  
  n <- length(x)
  
  if (n < 4L)
    stop("sample size must be at least 4")
  
  alternative <- match.arg(alternative)
  
  d  <- diff(x)
  mx <- mean(x)
  
  ss <- sum((x - mx)^2)
  if (ss <= 0)
    stop("'x' must contain at least two distinct observations")
  
  if (unbiased) {
    VN <- sum(d^2) / ss * n / (n - 1)
    Ex <- 2 * n / (n - 1)
    Vx <- 4 * n^2 * (n - 2) / ((n + 1) * (n - 1)^3)
  } else {
    VN <- sum(d^2) / ss
    Ex <- 2
    Vx <- 4 * (n - 2) / ((n + 1) * (n - 1))
  }
  
  z <- (VN - Ex) / sqrt(Vx)
  
  PVAL <- switch(
    alternative,
    "two.sided" = 2 * pnorm(-abs(z)),
    "less"      = pnorm(z),
    "greater"   = pnorm(z, lower.tail = FALSE)
  )
  
  structure(
    list(
      statistic   = c(z = z),
      parameter   = c(n = n),
      p.value     = as.numeric(PVAL),
      alternative = alternative,
      method      = "Von Neumann Successive Difference Test",
      data.name   = DNAME,
      vn          = VN
    ),
    class = "htest"
  )
}

