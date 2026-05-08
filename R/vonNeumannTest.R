
#' Von Neumann's Successive Difference Test 
#' 
#' A test for randomness or autocorrelation in a sequence, based on the mean 
#' square successive difference relative to the sample variance, closely 
#' related to the Durbin-Watson test.
#' 
#' A popular statistic to test for independence is the von Neumann ratio. 
#' 
#' The VN test statistic is in the unbiased case
#' \deqn{VN=\frac{\sum_{i=1}^{n-1}(x_i-x_{i+1})^2 \cdot
#' n}{\sum_{i=1}^{n}\left(x_i-\bar{x}\right)^2 \cdot (n-1)}
#' }{VN=\sum(x_i-x_{i+1})^2 / \sum(x_i-mean(x)^2 * n/n-1} It is known that
#' \eqn{(VN-\mu)/\sigma} is asymptotically standard normal, where
#' \eqn{\mu=\frac{2n}{n-1}}{\mu=2n/(n-1)} and \eqn{\sigma^2=4\cdot n^2
#' \frac{(n-2)}{(n+1)(n-1)^3}}{\sigma^2=[4*n^2 * (n-2)]/[(n+1)(n-1)^3]}.
#' 
#' The VN test statistic is in the original (biased) case
#' \deqn{VN=\frac{\sum_{i=1}^{n-1}(x_i-x_{i+1})^2}{\sum_{i=1}^{n}\left(x_i-\bar{x}\right)^2}}{VN=\sum(x_i-x_{i+1})^2
#' / \sum(x_i-mean(x)^2} The test statistic \eqn{(VN-2)/\sigma} is
#' asymptotically standard normal, where
#' \eqn{\sigma^2=\frac{4\cdot(n-2)}{(n+1)(n-1)}}{\sigma^2=[4*(n-2)]/[(n+1)(n-1)]}.
#' 
#' Missing values are silently removed.
#' 
#' @name vonNeumannTest
#' @param x a numeric vector containing the observations
#' @param alternative a character string specifying the alternative hypothesis,
#' must be one of \code{"two.sided"} (default), \code{"greater"} or
#' \code{"less"}. You can specify just the initial letter.
#' @param unbiased logical. In order for VN to be an unbiased estimate of the
#' true population value, the calculated value is multiplied by
#' \eqn{n/(n-1)}{n/(n-1)}. Default is TRUE. 
#' @return A list with class "htest" containing the components:
#' \item{statistic}{the value of the VN statistic and the normalized statistic
#' test.} \item{parameter, n}{the size of the data, after the remotion of
#' consecutive duplicate values.} \item{p.value}{the p-value of the test.}
#' \item{alternative}{a character string describing the alternative
#' hypothesis.} \item{method}{a character string indicating the test
#' performed.} \item{data.name}{a character string giving the name of the
#' data.}
#' 
#' @seealso \code{\link{bartelsRankTest}} 
#' @references von Neumann, J. (1941) Distribution of the ratio of the mean
#' square successive difference to the variance. \emph{Annals of Mathematical
#' Statistics} \bold{12}, 367-395.
#' 
#' @examples
#' 
#' set.seed(2)
#' vonNeumannTest(runif(20))


#' @rdname vonNeumannTest
#' @family test.correlation
#' @concept hypothesis-testing
#' @concept nonparametric
#' @concept time-series
#'
#'
#' @export
vonNeumannTest <- function (x, alternative = c("two.sided", "less", "greater"), unbiased=TRUE) {
  
  
  ## ToDo: use incomplete beta for exact p-values
  ## ************************
  ## see: von Neumann Successive Difference 1941
  ##
  # n <- 50
  # vx <- 1
  #
  # mu2 <- (4 * (3*n - 4)/(n-1)^2) * vx^2
  #
  # q2 <- (3*n^4 - 10*n^3 -18*n^2 + 79*n - 60) / (8*n^3 - 50*n + 48)
  # q1 <- (4 - mu2 * (q2 + 1) * (q2 + 3)) / (4 - mu2 * (q2 + 1))
  # a2 <- 2 * (q1 - q2 - 2) / (q2 + 1)
  # cc <- a2 ^(q1 - q2 - 2) / beta(q1 - q2 -1, q2+1)
  #
  # c(q1, q2, a2, cc)
  #
  # pbeta(0.75, shape1 = q1 - q2 -1, shape2= q2+1)
  # pbeta(0.75, shape1 = q1 - q2 -1, shape2= q2+1)
  #
  # beta(q1 - q2 -1, q2+1)
  
  
  alternative <- match.arg(alternative)
  
  dname <- deparse(substitute(x))
  
  x <- x[!is.na(x)]
  
  d <- diff(x)
  n <- length(x)
  mx <- mean(x)
  
  if(unbiased) {
    
    # http://www.chegg.com/homework-help/detecting-autocorrelation-von-neumann-ratio-test-assuming-re-chapter-12-problem-4-solution-9780073375779-exc
    
    VN <- sum(d^2) / sum((x - mx)^2) * n/(n-1)
    Ex <- 2 * n/(n-1)
    Vx <- 4 * n^2 * (n-2) / ((n+1) * (n-1)^3)
    z <- (VN - Ex) / sqrt(Vx)
    
  } else {
    VN <- sum(d^2) / sum((x - mx)^2)
    z <- (1-(VN/2)) / sqrt((n-2)/(n^2 - 1))
  }
  
  
  if (alternative == "less") {
    pval <- pnorm(z)
  }
  else if (alternative == "greater") {
    pval <- pnorm(z, lower.tail = FALSE)
  }
  else {
    pval <- 2 * pnorm(-abs(z))
  }
  names(VN) <- "VN"
  method <- "Von Neumann Successive Difference Test"
  
  rval <- list(statistic = c(VN, z=z), p.value = pval,
               method = method,
               alternative = alternative, data.name = dname,
               z = z)
  
  class(rval) <- "htest"
  return(rval)
  
}


