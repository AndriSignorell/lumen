
#' Bartels Rank Test of Randomness
#' 
#' A nonparametric test for randomness of a sequence (so it tests whether 
#' data is sampled randomly from an underlying population), based on the 
#' ratio of the mean square successive difference of ranks to the rank variance.
#' 
#' Data must at least be measured on an ordinal scale.
#' The RVN test statistic is
#' \deqn{RVN=\frac{\sum_{i=1}^{n-1}(R_i-R_{i+1})^2}{\sum_{i=1}^{n}\left(R_i-(n+1)/2\right)^2}}{RVN=\sum(R_i-R_{i+1})^2
#' / \sum(R_i-(n+1)/2)^2} where \eqn{R_i=rank(X_i), i=1,\dots,
#' n}{R_i=rank(X_i), i=1,...,n}. It is known that \eqn{(RVN-2)/\sigma} is
#' asymptotically standard normal, where
#' \eqn{\sigma^2=\frac{4(n-2)(5n^2-2n-9)}{5n(n+1)(n-1)^2}}{\sigma^2=[4(n-2)(5n^2-2n-9)]/[5n(n+1)(n-1)^2]}.
#' 
#' By using the alternative "\code{trend}" the null hypothesis of randomness is
#' tested against a trend. By using the alternative "\code{oscillation}" the
#' null hypothesis of randomness is tested against a systematic oscillation.
#' 
#' Missing values are silently removed.
#' 
#' Bartels test is a rank version of von Neumann's test.
#' 
#' @name bartelsRankTest
#' @param x a numeric vector containing the observations
#' @param alternative a character string specifying the alternative hypothesis,
#' must be one of "\code{two.sided}" (default), "\code{trend}" or
#' "\code{oscillation}".
#' @param method a character string specifying the method used to compute the
#' p-value. Must be one of \code{normal} (default), \code{beta} or \code{auto}.
#' @return A list with class "htest" containing the components:
#' \item{statistic}{the value of the normalized statistic test.}
#' \item{parameter, n}{the size of the data, after the remotion of consecutive
#' duplicate values.} \item{p.value}{the p-value of the test.}
#' \item{alternative}{a character string describing the alternative
#' hypothesis.} \item{method}{a character string indicating the test
#' performed.} \item{data.name}{a character string giving the name of the
#' data.} \item{rvn}{the value of the RVN statistic (not show on screen).}
#' \item{nm}{the value of the NM statistic, the numerator of RVN (not show on
#' screen).} \item{mu}{the mean value of the RVN statistic (not show on
#' screen).} \item{var}{the variance of the RVN statistic (not show on
#' screen).}
#' 
#' @note
#' Based on code by Frederico Caeiro.
#' 
#' @seealso \code{\link[randtests]{rank.test}}, \code{\link{runsTest}}
#' @references Bartels, R. (1982) The Rank Version of von Neumann's Ratio Test
#' for Randomness, \emph{Journal of the American Statistical Association},
#' \bold{77} (377), 40-46.
#' 
#' Gibbons, J.D. and Chakraborti, S. (2003) \emph{Nonparametric Statistical
#' Inference}, 4th ed. (pp. 97-98). URL:
#' \url{http://books.google.pt/books?id=dPhtioXwI9cC&lpg=PA97&ots=ZGaQCmuEUq}
#' 
#' von Neumann, J. (1941) Distribution of the ratio of the mean square
#' successive difference to the variance. \emph{Annals of Mathematical
#' Statistics} \bold{12}, 367-395.
#' 
#' @family topic.randomnessTests
#' @concept randomness
#' @concept rank-based
#' 
#' @examples
#' 
#' ## Example 5.1 in Gibbons and Chakraborti (2003), p.98.
#' ## Annual data on total number of tourists to the United States for 1970-1982.
#' 
#' years <- 1970:1982
#' tourists <- c(12362, 12739, 13057, 13955, 14123,  15698, 17523, 18610, 19842,
#'       20310, 22500, 23080, 21916)
#' plot(years, tourists, pch=20)
#' 
#' bartelsRankTest(tourists, alternative="trend", method="beta")
#' 
#' #  Bartels Ratio Test
#' #
#' # data:  tourists
#' # statistic = -3.6453, n = 13, p-value = 1.21e-08
#' # alternative hypothesis: trend
#' 
#' 
#' ## Example in Bartels (1982).
#' ## Changes in stock levels for 1968-1969 to 1977-1978 (in $A million), deflated by the
#' ## Australian gross domestic product (GDP) price index (base 1966-1967).
#' x <- c(528, 348, 264, -20, - 167, 575, 410, -4, 430, - 122)
#' 
#' bartelsRankTest(x, method="beta")



#' @rdname bartelsRankTest
#' @export
bartelsRankTest <- function(x, alternative = c("two.sided", "trend", "oscillation"),
                            method = c("normal", "beta", "auto")) {
  
  # Performs Bartels Ratio Test for Randomness.
  #
  # Args:
  #   x: a numeric vector containing the data.
  #   alternative hypothesis, must be one of "two.sided" (default), "left.sided" or "right.sided"
  #   pv.method: asymptotic aproximation method used to compute the p-value.
  #
  # Returns:
  #   statistic: the value of the RVN statistic test and the theoretical mean value and variance of the RVN statistic test.
  #   n: the sample size, after the remotion of consecutive duplicate values.
  #   p.value: the asymptotic p-value.
  #   method: a character string indicating the test performed.
  #   data.name: a character string giving the name of the data.
  #   alternative: a character string describing the alternative.
  
  pvalue <- match.arg(method, c("normal", "beta", "auto"))
  
  dname <- deparse(substitute(x))
  
  # Remove NAs
  x <- na.omit(x)
  stopifnot(is.numeric(x))
  n <- length(x)
  
  
  # if (alternative == "t"){alternative <- "two.sided"}
  # if (alternative == "l"){alternative <- "left.sided"}
  # if (alternative == "r"){alternative <- "right.sided"}
  # if (alternative != "two.sided" & alternative != "left.sided" & alternative != "right.sided")
  # {stop("must give a valid alternative")}
  
  alternative <- match.arg(alternative)
  
  if (n < 10){stop("sample size must be greater than 9")}
  
  # unique
  rk <- rank(x)
  d <- diff(rk)
  d.rank <- sum(rk^2) - n * (mean(rk)^2)
  RVN <- sum(d^2) / d.rank
  mu <- 2
  vr <- (4*(n-2)*(5*n^2-2*n-9))/(5*n*(n+1)*(n-1)^2)
  
  # Computes the p-value
  if (pvalue == "auto"){
    pvalue <- ifelse(n <= 100, "beta", "normal")
  }
  
  if (pvalue == "beta"){
    btp <- (5*n*(n+1)*(n-1)^2)/(2*(n-2)*(5*n^2-2*n-9))-1/2
    pv0 <- pbeta(RVN/4, shape1=btp, shape2=btp)
  }
  if (pvalue=="normal"){
    pv0 <- pnorm((RVN - mu) / sqrt(vr))
  }
  
  if (alternative=="two.sided"){
    pv <- 2 * min(pv0, 1 - pv0)
    alternative <- "nonrandomness"
  }
  if (alternative == "trend"){
    pv <- pv0
    alternative <- "trend"
  }
  if (alternative == "oscillation"){
    pv <- 1 - pv0
    alternative <- "systematic oscillation"
  }
  
  test <- (RVN - mu) / sqrt(vr)
  rval <- list(statistic = c(RVN=RVN, z=test), nm=sum(d^2), rvn=RVN, mu=mu, var=vr, p.value = pv,
               method = "Bartels Ratio Test", data.name = dname, parameter=c(n=n), n=n, alternative=alternative)
  class(rval) <- "htest"
  return(rval)
  
}

