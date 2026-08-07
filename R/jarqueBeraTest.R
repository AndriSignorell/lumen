
#' (Robust) Jarque-Bera Test for Assessing Normality From Skewness and Kurtosis
#'
#' A goodness-of-fit test for normality based on sample skewness and
#' kurtosis, in its classical form or the robust modification of Gel and
#' Gastwirth (2008).
#'
#' This function performs either the classical Jarque-Bera test or the
#' robust Jarque-Bera test proposed by Gel and Gastwirth (2008). The robust
#' version (the default) replaces the classical standard deviation by the
#' average absolute deviation from the median (MAAD), which makes the test
#' considerably less sensitive to outliers. The kurtosis component is then
#' scaled with the constant \eqn{C_2 = 64} derived in Gel and Gastwirth
#' (2008), the skewness component with \eqn{C_1 = 6} as in the classical
#' test.
#'
#' With \code{method = "chisq"} the statistic is referred to its asymptotic
#' chi-squared distribution with 2 degrees of freedom. With
#' \code{method = "mc"} the p-value is estimated by Monte Carlo simulation
#' from \code{R} samples of the same size drawn from the standard normal
#' distribution, using the finite-sample correction
#' \eqn{(m + 1)/(R + 1)}, where \eqn{m} is the number of simulated
#' statistics at least as large as the observed one.
#'
#' Missing values are silently removed.
#'
#' @param x a numeric vector of data values.
#' @param robust logical, if \code{TRUE} (default) the robust test of Gel
#' and Gastwirth (2008) is performed, otherwise the classical Jarque-Bera
#' test.
#' @param method a character string specifying how the p-value is computed,
#' one of \code{"chisq"} (default, asymptotic chi-squared approximation) or
#' \code{"mc"} (Monte Carlo simulation).
#' @param R the number of Monte Carlo replicates used for
#' \code{method = "mc"} (default is \code{1000}).
#' @return A list with class \code{"htest"} containing the following
#' components:
#' \item{statistic}{the value of the test statistic.}
#' \item{parameter}{the degrees of freedom (\code{method = "chisq"}), or
#' the number of Monte Carlo replicates (\code{method = "mc"}).}
#' \item{p.value}{the p-value of the test.}
#' \item{method}{a character string indicating the test performed and how
#' the p-value was computed.}
#' \item{data.name}{a character string giving the name of the data.}
#'
#' @note
#' Based on code by W. Wallace Hui, Yulia R. Gel, Joseph L. Gastwirth and
#' Weiwen Miao previously published as \code{rjb.test()} in the
#' \pkg{lawstat} package, adapted to conform to package standards.
#'
#' @references
#' Gel, Y. R. and Gastwirth, J. L. (2008) A robust modification of the
#' Jarque-Bera test of normality. \emph{Economics Letters}, 99, 30-32.
#'
#' Jarque, C. and Bera, A. (1980) Efficient tests for normality,
#' homoscedasticity and serial independence of regression residuals.
#' \emph{Economics Letters}, 6, 255-259.
#'
#' @seealso [shapiro.test()]
#'
#' @examples
#' set.seed(1)
#' x <- rnorm(100)
#'
#' jarqueBeraTest(x)                   # robust version (default)
#' jarqueBeraTest(x, robust = FALSE)   # classical Jarque-Bera test
#'
#' # Monte Carlo p-value
#' jarqueBeraTest(x, method = "mc", R = 5000)
#'
#' # heavy-tailed alternative
#' jarqueBeraTest(rt(100, df = 2))
#'
#' @family test.normality
#' @concept normality-test
#' @concept goodness-of-fit
#'
#' @export
jarqueBeraTest <- function(x,
                           robust = TRUE,
                           method = c("chisq", "mc"),
                           R = 1000) {

  DNAME <- deparse1(substitute(x))

  method <- match.arg(method)

  if (!is.null(dim(x)))
    stop("'x' must be a numeric vector or univariate time series")

  x <- as.numeric(x)
  x <- x[!is.na(x)]
  n <- length(x)

  if (n < 3L)
    stop("sample size must be at least 3")

  if (method == "mc") {
    R <- as.integer(R)
    if (is.na(R) || R < 1L)
      stop("'R' must be a positive integer")
  }

  STATISTIC <- .jbStatistic(x, robust = robust)

  PVAL <- switch(
    method,
    chisq = pchisq(STATISTIC, df = 2, lower.tail = FALSE),
    mc    = .jbPvalueMC(statistic = STATISTIC, n = n,
                        robust = robust, R = R)
  )

  METHOD <- paste0(
    if (robust) "Robust " else "",
    "Jarque-Bera test (",
    if (method == "mc") "Monte Carlo" else "chi-squared approximation",
    ")"
  )

  PARAMETER <- switch(method,
                      chisq = c(df = 2),
                      mc    = c(R = R))

  structure(
    list(
      statistic = c("X-squared" = STATISTIC),
      parameter = PARAMETER,
      p.value   = PVAL,
      method    = METHOD,
      data.name = DNAME
    ),
    class = "htest"
  )
}


# == internal helper functions ============================================


.jbStatistic <- function(x, robust = TRUE) {

  if (diff(range(x)) == 0)
    stop("all values are identical")

  n <- length(x)

  zc <- x - mean(x)

  m2 <- mean(zc^2)
  m3 <- mean(zc^3)
  m4 <- mean(zc^4)

  if (robust) {

    # normalising constant such that
    # E[sqrt(pi / 2) * |Z - median(Z)|] = sigma for Z ~ N(0, sigma^2)
    J <- sqrt(pi / 2) * mean(abs(x - median(x)))

    if (J <= 0)
      stop("robust scale estimate is zero")

    sq_skew  <- (m3 / J^3)^2
    kurtosis <- m4 / J^4

    # Gel and Gastwirth (2008): C1 = 6, C2 = 64
    vs <- 6 / n
    vk <- 64 / n

  } else {

    sq_skew  <- (m3 / m2^(3 / 2))^2
    kurtosis <- m4 / m2^2

    vs <- 6 / n
    vk <- 24 / n
  }

  sq_skew / vs + (kurtosis - 3)^2 / vk
}


.jbPvalueMC <- function(statistic, n, robust, R) {

  sim <- replicate(R, .jbStatistic(rnorm(n), robust = robust))

  # Monte Carlo p-value with finite-sample correction
  (sum(sim >= statistic) + 1) / (R + 1)
}
