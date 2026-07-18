
#' Augmented Dickey-Fuller Unit Root Test
#'
#' Performs the augmented Dickey-Fuller test for a unit root in a time
#' series, testing the null hypothesis of nonstationarity against the
#' alternative of (trend-)stationarity. It extends the simple Dickey-Fuller
#' test by including lagged difference terms to account for autocorrelation.
#'
#' If \code{type} is set to \code{"none"} neither an intercept nor a trend
#' is included in the test regression, \code{"drift"} adds an intercept and
#' \code{"trend"} adds both an intercept and a linear trend.
#'
#' The reported test statistic is the t statistic of the lagged level
#' (\code{tau1}, \code{tau2} or \code{tau3}, depending on \code{type}).
#' For \code{type = "drift"} and \code{type = "trend"} the F statistics
#' \code{phi1} resp. \code{phi2} and \code{phi3} of Dickey and Fuller (1981),
#' testing joint hypotheses on the deterministic terms, are reported as
#' additional statistics; their critical values are contained in the
#' \code{critical.values} component of the result.
#'
#' The p-value refers to the tau statistic and is obtained by linear
#' interpolation in the finite sample quantiles given by Fuller (1976),
#' following the approach of \code{tseries::adf.test()}. If the statistic
#' falls outside the range of the table, the p-value is reported as the
#' respective boundary (0.01 or 0.99) and a warning is issued. The critical
#' values are taken from Hamilton (1994) and Dickey and Fuller (1981).
#'
#' Missing values are not allowed.
#'
#' @param y numeric vector or univariate time series to be tested for a
#' unit root.
#' @param type the deterministic part of the test regression, one of
#' \code{"none"} (default), \code{"drift"} or \code{"trend"}.
#' @param lags the number of lagged difference terms to be included. If
#' \code{selectLags} is not \code{"fixed"}, this is the maximum number of
#' lags considered in the lag selection.
#' @param selectLags the lag selection method, one of \code{"fixed"}
#' (default, use \code{lags} as given), \code{"aic"} or \code{"bic"}
#' (choose the lag order up to \code{lags} that minimizes the respective
#' information criterion). Case-insensitive, so \code{"AIC"} and
#' \code{"BIC"} are accepted as well.
#' 
#' @return An object of class \code{"htest"} containing the following
#' components:
#'   \item{\code{statistic}}{the tau statistic, followed by the phi
#'     statistic(s) if \code{type} is \code{"drift"} or \code{"trend"}. The
#'     p-value refers to the tau statistic.}
#'   \item{\code{parameter}}{the number of lagged differences included in the
#'     test regression (after lag selection, if requested).}
#'   \item{\code{p.value}}{the interpolated p-value of the tau statistic.}
#'   \item{\code{critical.values}}{a matrix with the 1\%, 5\% and 10\% critical
#'     values of all reported statistics, interpolated for the effective
#'     sample size (not shown on screen).}
#'   \item{\code{alternative}}{a character string describing the alternative
#'     hypothesis.}
#'   \item{\code{method}}{a character string indicating the test performed.}
#'   \item{\code{data.name}}{a character string giving the name of the data.}
#'
#' @note
#' Based on code by Bernhard Pfaff previously published in the \pkg{urca}
#' package, adapted to conform to package standards.
#'
#' @references Dickey, D. A. and Fuller, W. A. (1979) Distribution of the
#' estimators for autoregressive time series with a unit root, \emph{Journal
#' of the American Statistical Association}, \bold{74}, 427--431.
#'
#' Dickey, D. A. and Fuller, W. A. (1981) Likelihood ratio statistics for
#' autoregressive time series with a unit root, \emph{Econometrica},
#' \bold{49}, 1057--1072.
#'
#' Fuller, W. A. (1976) \emph{Introduction to Statistical Time Series},
#' New York: Wiley.
#'
#' Hamilton, J. D. (1994) \emph{Time Series Analysis}, Princeton:
#' Princeton University Press.
#'
#' @seealso \code{\link{kpssTest}}
#'
#' @examples
#' adfTest(AirPassengers, lags = 3, type = "trend")
#'
#' # a random walk should not be rejected
#' set.seed(5)
#' rw <- cumsum(rnorm(200))
#' adfTest(rw, type = "drift")
#'
#' @family test.timeseries
#' @concept unit-root
#' @concept stationarity
#'
#' @export
adfTest <- function(y, type = c("none", "drift", "trend"),
                    lags = 1, selectLags = c("fixed", "aic", "bic")) {

  type <- match.arg(type)
  selectLags <- match.arg(tolower(selectLags), c("fixed", "aic", "bic"))

  DNAME <- deparse(substitute(y))

  if (ncol(as.matrix(y)) > 1)
    stop("'y' must be a vector or a univariate time series")
  if (anyNA(y))
    stop("NAs in 'y'")

  lag <- as.integer(lags)
  if (is.na(lag) || lag < 0)
    stop("'lags' must be a non-negative integer")

  y <- as.vector(y)
  lags <- lag + 1L

  z <- diff(y)
  n <- length(z)
  if (n < lags + 2L)
    stop("'y' is too short for the requested number of lags")

  x       <- embed(z, lags)
  z.diff  <- x[, 1]
  z.lag.1 <- y[lags:n]
  tt      <- lags:n

  # test regression with k - 1 lagged differences (k = embed order)
  fitADF <- function(k) {
    if (k > 1) {
      z.diff.lag <- x[, 2:k]
      switch(type,
             none  = lm(z.diff ~ z.lag.1 - 1 + z.diff.lag),
             drift = lm(z.diff ~ z.lag.1 + 1 + z.diff.lag),
             trend = lm(z.diff ~ z.lag.1 + 1 + tt + z.diff.lag))
    } else {
      switch(type,
             none  = lm(z.diff ~ z.lag.1 - 1),
             drift = lm(z.diff ~ z.lag.1 + 1),
             trend = lm(z.diff ~ z.lag.1 + 1 + tt))
    }
  }

  # lag selection by information criterion
  if (lags > 1 && selectLags != "fixed") {
    penalty <- switch(selectLags, aic = 2, bic = log(length(z.diff)))
    critRes <- rep(NA_real_, lags)
    for (i in 2:lags)
      critRes[i] <- AIC(fitADF(i), k = penalty)
    lags <- which.min(critRes)
  }

  result <- fitADF(lags)
  tau <- coef(summary(result))[if (type == "none") 1 else 2, 3]

  # phi F statistics (Dickey and Fuller, 1981) for the joint hypotheses
  # on the deterministic terms
  phi <- NULL
  if (type != "none") {
    z.diff.lag <- if (lags > 1) x[, 2:lags] else NULL
    if (type == "drift") {
      phi1.reg <- if (lags > 1) lm(z.diff ~ -1 + z.diff.lag) else
        lm(z.diff ~ -1)
      phi <- c(phi1 = anova(phi1.reg, result)$F[2])
    } else {
      phi2.reg <- if (lags > 1) lm(z.diff ~ -1 + z.diff.lag) else
        lm(z.diff ~ -1)
      phi3.reg <- if (lags > 1) lm(z.diff ~ z.diff.lag) else
        lm(z.diff ~ 1)
      phi <- c(phi2 = anova(phi2.reg, result)$F[2],
               phi3 = anova(phi3.reg, result)$F[2])
    }
  }

  tauname <- switch(type, none = "tau1", drift = "tau2", trend = "tau3")
  STATISTIC <- c(tau, phi)
  names(STATISTIC)[1] <- tauname

  # p-value for tau by interpolation in Fuller's table,
  # following tseries::adf.test()
  tauQuant <- apply(.adfTauTable[[type]], 2,
                    function(cv) approx(.adfTabN, cv, n, rule = 2)$y)
  PVAL <- approx(tauQuant, .adfTauProbs, tau, rule = 2)$y
  if (is.na(approx(tauQuant, .adfTauProbs, tau, rule = 1)$y)) {
    if (PVAL == min(.adfTauProbs))
      warning("p-value smaller than reported p-value")
    else
      warning("p-value greater than reported p-value")
  }

  # critical values of all reported statistics at the effective sample size
  cvals <- rbind(tauQuant[c(1, 3, 4)],
                 if (type != "none")
                   t(sapply(.adfPhiTable[names(phi)], function(tab)
                     apply(tab, 2, function(cv)
                       approx(.adfTabN, cv, n, rule = 2)$y))))
  dimnames(cvals) <- list(names(STATISTIC), c("1pct", "5pct", "10pct"))

  structure(
    list(statistic       = STATISTIC,
         parameter       = c(lags = lags - 1L),
         p.value         = PVAL,
         critical.values = cvals,
         alternative     = switch(type, trend = "trend-stationary",
                                  "stationary"),
         method          = "Augmented Dickey-Fuller Test",
         data.name       = DNAME),
    class = "htest")
}



# == internal constants ================================================


# Empirical quantiles of the Dickey-Fuller tau statistics,
# Fuller (1976), Table 8.5.2 (reproduced in Hamilton 1994, Table B.6).
# Rows: sample sizes 25, 50, 100, 250, 500, Inf;
# columns: probabilities 0.01, 0.025, 0.05, 0.10, 0.90, 0.95, 0.975, 0.99.

.adfTauProbs <- c(0.01, 0.025, 0.05, 0.10, 0.90, 0.95, 0.975, 0.99)
.adfTabN     <- c(25, 50, 100, 250, 500, 1e5)

.adfTauTable <- list(
  none = rbind(
    c(-2.66, -2.26, -1.95, -1.60, 0.92, 1.33, 1.70, 2.16),
    c(-2.62, -2.25, -1.95, -1.61, 0.91, 1.31, 1.66, 2.08),
    c(-2.60, -2.24, -1.95, -1.61, 0.90, 1.29, 1.64, 2.03),
    c(-2.58, -2.23, -1.95, -1.62, 0.89, 1.29, 1.63, 2.01),
    c(-2.58, -2.23, -1.95, -1.62, 0.89, 1.28, 1.62, 2.00),
    c(-2.58, -2.23, -1.95, -1.62, 0.89, 1.28, 1.62, 2.00)),
  drift = rbind(
    c(-3.75, -3.33, -3.00, -2.63, -0.37,  0.00, 0.34, 0.72),
    c(-3.58, -3.22, -2.93, -2.60, -0.40, -0.03, 0.29, 0.66),
    c(-3.51, -3.17, -2.89, -2.58, -0.42, -0.05, 0.26, 0.63),
    c(-3.46, -3.14, -2.88, -2.57, -0.42, -0.06, 0.24, 0.62),
    c(-3.44, -3.13, -2.87, -2.57, -0.43, -0.07, 0.24, 0.61),
    c(-3.43, -3.12, -2.86, -2.57, -0.44, -0.07, 0.23, 0.60)),
  trend = rbind(
    c(-4.38, -3.95, -3.60, -3.24, -1.14, -0.80, -0.50, -0.15),
    c(-4.15, -3.80, -3.50, -3.18, -1.19, -0.87, -0.58, -0.24),
    c(-4.04, -3.73, -3.45, -3.15, -1.22, -0.90, -0.62, -0.28),
    c(-3.99, -3.69, -3.43, -3.13, -1.23, -0.92, -0.64, -0.31),
    c(-3.98, -3.68, -3.42, -3.13, -1.24, -0.93, -0.65, -0.32),
    c(-3.96, -3.66, -3.41, -3.12, -1.25, -0.94, -0.66, -0.33)))

# Critical values (1%, 5%, 10%) of the phi F statistics,
# Dickey and Fuller (1981), Tables IV-VI; rows as above.
# All values verified against the original paper. Note that urca
# carries a transcription error in the phi3 row for n = 250
# (5%/10% values duplicated from the n = 100 row); the correct
# values per Table VI are 6.34 and 5.39.

.adfPhiTable <- list(
  phi1 = rbind(
    c(7.88, 5.18, 4.12),
    c(7.06, 4.86, 3.94),
    c(6.70, 4.71, 3.86),
    c(6.52, 4.63, 3.81),
    c(6.47, 4.61, 3.79),
    c(6.43, 4.59, 3.78)),
  phi2 = rbind(
    c(8.21, 5.68, 4.67),
    c(7.02, 5.13, 4.31),
    c(6.50, 4.88, 4.16),
    c(6.22, 4.75, 4.07),
    c(6.15, 4.71, 4.05),
    c(6.09, 4.68, 4.03)),
  phi3 = rbind(
    c(10.61, 7.24, 5.91),
    c( 9.31, 6.73, 5.61),
    c( 8.73, 6.49, 5.47),
    c( 8.43, 6.34, 5.39),
    c( 8.34, 6.30, 5.36),
    c( 8.27, 6.25, 5.34)))

