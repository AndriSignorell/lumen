
#' Kwiatkowski-Phillips-Schmidt-Shin Test for Stationarity
#'
#' A test for stationarity in time series, complementary to unit root tests
#' such as the ADF test: it tests the null hypothesis of stationarity
#' against the alternative of a unit root.
#'
#' Performs the KPSS test (Kwiatkowski et al., 1992), where the null
#' hypothesis is stationarity. The test types specify as deterministic
#' component either a constant (\code{"mu"}, null hypothesis of level
#' stationarity) or a constant with linear trend (\code{"tau"}, null
#' hypothesis of trend stationarity).
#'
#' \code{lags = "short"} sets the number of lags used for the long-run
#' variance estimation to \eqn{4 (n/100)^{1/4}}{4*(n/100)^(1/4)}, whereas
#' \code{lags = "long"} sets it to \eqn{12 (n/100)^{1/4}}{12*(n/100)^(1/4)}
#' (each truncated to an integer). With \code{lags = "nil"} no error
#' correction is made. Alternatively, an explicit number of lags can be
#' given via \code{useLag}, which then takes precedence.
#'
#' The p-value is obtained by linear interpolation in the asymptotic
#' critical values of Kwiatkowski et al. (1992, Table 1), following the
#' approach of \code{tseries::kpss.test()}. If the statistic falls outside
#' the range of the table, the p-value is reported as the respective
#' boundary (0.01 or 0.10) and a warning is issued.
#'
#' Missing values are silently removed.
#'
#' @param y numeric vector or univariate time series to be tested for
#' stationarity.
#' @param type the deterministic part of the model, one of \code{"mu"}
#' (default, constant) or \code{"tau"} (constant plus linear trend).
#' @param lags the rule for the number of lags used for the error term
#' correction, one of \code{"short"} (default), \code{"long"} or
#' \code{"nil"}. See the Details. Ignored if \code{useLag} is given.
#' @param useLag an optional integer explicitly specifying the number of
#' lags, overriding \code{lags}.
#'
#' @return An object of class \code{"htest"} containing the following
#' components:
#' \item{statistic}{the value of the KPSS test statistic.}
#' \item{parameter}{the number of lags used for the error term correction.}
#' \item{p.value}{the interpolated p-value of the test.}
#' \item{critical.values}{the asymptotic critical values at the 10%, 5%,
#' 2.5% and 1% significance levels (not shown on screen).}
#' \item{alternative}{a character string describing the alternative
#' hypothesis.}
#' \item{method}{a character string indicating the test performed.}
#' \item{data.name}{a character string giving the name of the data.}
#'
#' @note
#' Based on code by Bernhard Pfaff previously published in the \pkg{urca}
#' package, adapted to conform to package standards.
#'
#' @references Kwiatkowski, D., Phillips, P. C. B., Schmidt, P. and
#' Shin, Y. (1992) Testing the null hypothesis of stationarity against the
#' alternative of a unit root: How sure are we that economic time series
#' have a unit root? \emph{Journal of Econometrics}, \bold{54}, 159--178.
#'
#' @seealso \code{\link{adfTest}}
#'
#' @examples
#' # trend-stationary series: null hypothesis is not rejected
#' set.seed(1)
#' x <- 0.2 * seq_len(200) + rnorm(200)
#' kpssTest(x, type = "tau")
#'
#' # random walk: null hypothesis of level stationarity is rejected
#' set.seed(2)
#' rw <- cumsum(rnorm(200))
#' kpssTest(rw, type = "mu")
#'
#' kpssTest(AirPassengers, type = "tau", lags = "short")
#'
#' @family test.timeseries
#' @concept stationarity
#' @concept unit-root
#'
#' @export
kpssTest <- function(y, type = c("mu", "tau"),
                     lags = c("short", "long", "nil"), useLag = NULL) {

  DNAME <- deparse1(substitute(y))

  type <- match.arg(type)
  lags <- match.arg(lags)

  y <- as.vector(y)
  y <- y[!is.na(y)]
  n <- length(y)

  # number of lags for the error term correction
  lmax <- switch(lags,
                 short = trunc(4 * (n / 100)^0.25),
                 long  = trunc(12 * (n / 100)^0.25),
                 nil   = 0L)

  if (!is.null(useLag)) {
    useLag <- suppressWarnings(as.integer(useLag))
    if (length(useLag) != 1L || is.na(useLag) || useLag < 0L)
      warning("'useLag' must be a single non-negative integer; ",
              "using lags = '", lags, "' instead")
    else
      lmax <- useLag
  }

  # asymptotic critical values, Kwiatkowski et al. (1992), Table 1
  probs <- c(0.10, 0.05, 0.025, 0.01)
  cval  <- switch(type,
                  mu  = c(0.347, 0.463, 0.574, 0.739),
                  tau = c(0.119, 0.146, 0.176, 0.216))
  cvals <- matrix(cval, nrow = 1,
                  dimnames = list("critical values",
                                  c("10pct", "5pct", "2.5pct", "1pct")))

  # residuals of the deterministic part
  res <- switch(type,
                mu  = y - mean(y),
                tau = residuals(lm(y ~ seq_len(n))))

  # KPSS statistic: partial sums over long-run variance estimate
  S  <- cumsum(res)
  s2 <- sum(res^2) / n

  denominator <- if (lmax == 0) {
    s2
  } else {
    index    <- seq_len(lmax)
    x.cov    <- vapply(index, function(k)
      sum(res[-seq_len(k)] * res[seq_len(n - k)]), numeric(1))
    bartlett <- 1 - index / (lmax + 1)
    s2 + 2 / n * sum(bartlett * x.cov)
  }

  STATISTIC <- (sum(S^2) / n^2) / denominator

  # p-value by interpolation in the critical values,
  # following tseries::kpss.test()
  PVAL <- approx(cval, probs, STATISTIC, rule = 2)$y
  if (is.na(approx(cval, probs, STATISTIC, rule = 1)$y)) {
    if (PVAL == max(probs))
      warning("p-value greater than reported p-value")
    else
      warning("p-value smaller than reported p-value")
  }

  structure(
    list(
      statistic       = c(KPSS = as.numeric(STATISTIC)),
      parameter       = c(lags = as.integer(lmax)),
      p.value         = PVAL,
      critical.values = cvals,
      alternative     = switch(type,
                               mu  = "not level stationary (unit root)",
                               tau = "not trend stationary (unit root)"),
      method          = paste("KPSS test for",
                              switch(type, mu = "level", tau = "trend"),
                              "stationarity"),
      data.name       = DNAME
    ),
    class = "htest"
  )
}
