
#' Durbin-Watson Test
#'
#' Tests for first-order autocorrelation in the residuals of a linear
#' regression model, based on the ratio of successive squared residual
#' differences to the total residual sum of squares.
#'
#' @description
#' The Durbin-Watson test has the null hypothesis that the autocorrelation of
#' the disturbances is 0. The alternative hypothesis can be specified as
#' greater than, not equal to, or less than 0 via the \code{alternative}
#' argument.
#'
#' Under the assumption of normally distributed disturbances, the null
#' distribution of the Durbin-Watson statistic is the distribution of a linear
#' combination of chi-squared variables. The p-value is computed using the
#' \code{"pan"} algorithm (Farebrother, 1980, 1984), implemented via Rcpp.
#' For large sample sizes the algorithm might fail to compute the p-value;
#' in that case a warning is printed and an approximate p-value is returned
#' instead, computed via a normal approximation using the mean and variance of
#' the Durbin-Watson statistic.
#'
#' Three methods are dispatched:
#' \describe{
#'   \item{\code{formula}}{Fits the model from scratch using \code{data}.}
#'   \item{\code{lm}}{Extracts design matrix and response from a fitted
#'     \code{"lm"} object. The \code{data.name} field in the result reflects
#'     the model formula, not the name of the object.}
#'   \item{\code{numeric}}{Treats \code{x} as a vector of pre-computed
#'     residuals and fits an intercept-only design matrix. Note that this is
#'     \emph{not} equivalent to testing residuals from a fitted model with
#'     predictors; use the \code{lm} or \code{formula} method when a model
#'     exists.}
#' }
#'
#' @param x a symbolic description of the model to be tested (a
#'   \code{formula}), a fitted \code{"lm"} object, or a numeric vector of
#'   residuals.
#' @param data an optional data frame containing the variables in the model.
#'   Only used for the \code{formula} method. By default variables are taken
#'   from the environment from which \code{durbinWatsonTest} is called.
#' @param orderBy either a numeric vector or a one-sided formula
#'   \code{~ z}. Observations are reordered by increasing \code{z} before
#'   the test is applied. If \code{NULL} (the default) observations are
#'   assumed to be already ordered (e.g. a time series).
#' @param alternative a character string specifying the alternative hypothesis,
#'   must be one of \code{"greater"} (default), \code{"two.sided"}, or
#'   \code{"less"}.
#' @param iterations an integer specifying the number of iterations used by
#'   the \code{"pan"} algorithm when computing the exact p-value.
#' @param exact logical. If \code{TRUE} the exact p-value is computed via the
#'   \code{"pan"} algorithm; if \code{FALSE} a normal approximation is used.
#'   The default is \code{TRUE} for sample sizes below 100 and \code{FALSE}
#'   otherwise.
#' @param tol numeric tolerance. Eigenvalues smaller than \code{tol} are
#'   treated as zero.
#' @param ... further arguments passed to or from other methods.
#'
#' @return An object of class \code{"htest"} with the following components:
#' \describe{
#'   \item{\code{statistic}}{the Durbin-Watson test statistic.}
#'   \item{\code{p.value}}{the corresponding p-value.}
#'   \item{\code{alternative}}{a character string describing the alternative
#'     hypothesis.}
#'   \item{\code{method}}{a character string identifying the test.}
#'   \item{\code{data.name}}{a character string describing the data.}
#' }
#'
#' @seealso \code{\link{lm}}
#'
#' @author Torsten Hothorn, Achim Zeileis, Clint Cummins (original
#'   \code{"pan"} algorithm); Andri Signorell (Rcpp reimplementation);
#'   Giovanni Millo, David Mitchell (further contributions)
#'
#' @references
#' J. Durbin & G.S. Watson (1950), Testing for Serial Correlation in Least
#' Squares Regression I. \emph{Biometrika} \bold{37}, 409--428.
#'
#' J. Durbin & G.S. Watson (1951), Testing for Serial Correlation in Least
#' Squares Regression II. \emph{Biometrika} \bold{38}, 159--178.
#'
#' J. Durbin & G.S. Watson (1971), Testing for Serial Correlation in Least
#' Squares Regression III. \emph{Biometrika} \bold{58}, 1--19.
#'
#' R.W. Farebrother (1980), Pan's Procedure for the Tail Probabilities of the
#' Durbin-Watson Statistic. \emph{Applied Statistics} \bold{29}, 224--227.
#' (C++ reimplementation by Andri Signorell.)
#'
#' R.W. Farebrother (1984), AS R53 A Remark on Algorithms AS 106, AS 153 and
#' AS 155. \emph{Applied Statistics} \bold{33}, 366--369.
#'
#' W. Krämer & H. Sonnberger (1986), \emph{The Linear Regression Model under
#' Test}. Heidelberg: Physica.
#'
#' J. Racine & R. Hyndman (2002), Using R To Teach Econometrics.
#' \emph{Journal of Applied Econometrics} \bold{17}, 175--189.
#' 
#' 
#' @examples
#' ## formula method
#' x <- rep(c(-1, 1), 50)
#'
#' err1 <- rnorm(100)
#' durbinWatsonTest(y ~ x, data = data.frame(y = 1 + x + err1, x = x))
#'
#' ## autocorrelated errors (rho = 0.9)
#' err2 <- stats::filter(err1, 0.9, method = "recursive")
#' durbinWatsonTest(y ~ x, data = data.frame(y = 1 + x + err2, x = x))
#'
#' ## lm method
#' fit <- lm(y ~ x, data = data.frame(y = 1 + x + err1, x = x))
#' durbinWatsonTest(fit)
#'
#' ## numeric method (pre-computed residuals — intercept-only design assumed)
#' e_t <- c(-32.33, -26.603, 2.215, -16.967, -1.148, -2.512, -1.967, 11.669,
#'          -0.513, 27.032, -4.422, 40.032, 23.577, 33.94, -2.787, -8.606,
#'           0.575, 6.848, -18.971, -29.063)
#' durbinWatsonTest(e_t)
#'
#' @family test.correlation
#' @concept hypothesis-testing
#' @concept regression
#' @concept time-series


# ── Generic ───────────────────────────────────────────────────────────────────

#' @rdname durbinWatsonTest
#' @export
durbinWatsonTest <- function(x, ...) UseMethod("durbinWatsonTest")


# ── Methoden ──────────────────────────────────────────────────────────────────

#' @rdname durbinWatsonTest
#' @export
durbinWatsonTest.formula <- function(x, data = list(), orderBy = NULL,
                                     alternative = c("greater", "two.sided", "less"),
                                     iterations = 15, exact = NULL,
                                     tol = 1e-10, ...) {
  mf  <- model.frame(x, data = data)
  y   <- model.response(mf)
  X   <- model.matrix(attr(mf, "terms"), mf)          # (3) gleicher model frame
  
  .dw_order(
    X           = X,
    y           = y,
    orderBy     = orderBy,
    data        = data,
    dname       = paste(deparse(substitute(x))),
    alternative = alternative,
    iterations  = iterations,
    exact       = exact,
    tol         = tol
  )
}

#' @rdname durbinWatsonTest
#' @export
durbinWatsonTest.lm <- function(x, orderBy = NULL,
                                alternative = c("greater", "two.sided", "less"),
                                iterations = 15, exact = NULL,
                                tol = 1e-10, ...) {
  if (!is.null(w <- weights(x)))
    if (!isTRUE(all.equal(as.vector(w), rep(1L, length(w)))))
      stop("weighted regressions are not supported")
  
  X <- if (is.matrix(x$x)) x$x else model.matrix(terms(x), model.frame(x))
  y <- if (is.vector(x$y)) x$y else model.response(model.frame(x))
  
  .dw_order(
    X           = X,
    y           = y,
    orderBy     = orderBy,
    data        = list(),
    dname       = paste(deparse(formula(x)), collapse = ""), # (2) Formel, nicht Objektname
    alternative = alternative,
    iterations  = iterations,
    exact       = exact,
    tol         = tol
  )
}

#' @rdname durbinWatsonTest
#' @export
durbinWatsonTest.numeric <- function(x, orderBy = NULL,                  # (1) explizit numeric
                                     alternative = c("greater", "two.sided", "less"),
                                     iterations = 15, exact = NULL,
                                     tol = 1e-10, ...) {
  X <- matrix(1, nrow = length(x), ncol = 1)
  
  .dw_order(
    X           = X,
    y           = x,
    orderBy     = orderBy,
    data        = list(),
    dname       = paste(deparse(substitute(x))),
    alternative = alternative,
    iterations  = iterations,
    exact       = exact,
    tol         = tol
  )
}

#' @rdname durbinWatsonTest
#' @export
durbinWatsonTest.default <- function(x, ...) {                           # (1) sprechender Fehler
  stop("no applicable method for objects of class ", sQuote(class(x)[1]))
}


# ── Interne Hilfsfunktionen ───────────────────────────────────────────────────

.dw_order <- function(X, y, orderBy, data, dname,
                      alternative, iterations, exact, tol) {
  if (!is.null(orderBy)) {
    mm <- model.matrix(orderBy, data = data)           # (4) kein doppelter Aufruf
    z  <- mm[, ncol(mm)]
    X  <- as.matrix(X[order(z), ])
    y  <- y[order(z)]
  }
  .dw_compute(
    X           = X,
    y           = y,
    dname       = dname,
    alternative = alternative,
    iterations  = iterations,
    exact       = exact,
    tol         = tol
  )
}

.dw_compute <- function(X, y, dname,
                        alternative = c("greater", "two.sided", "less"),
                        iterations = 15, exact = NULL, tol = 1e-10) {
  alternative <- match.arg(alternative)
  n           <- nrow(X)
  k           <- ncol(X)
  if (is.null(exact)) exact <- (n < 100)
  
  res     <- lm.fit(X, y)$residuals
  dw      <- sum(diff(res)^2) / sum(res^2)
  XtX_inv <- chol2inv(qr.R(qr(X)))                    # (6) sprechender Name
  
  pval <- if (n < 3) {
    warning("not enough observations for computing a p value, set to 1")
    1
  } else {
    .dw_pvalue(
      dw          = dw,
      X           = X,
      XtX_inv     = XtX_inv,
      n           = n,
      k           = k,
      alternative = alternative,
      iterations  = iterations,
      exact       = exact,
      tol         = tol
    )
  }
  
  alternative_str <- switch(alternative,
                            "two.sided" = "true autocorrelation is not 0",
                            "less"      = "true autocorrelation is less than 0",
                            "greater"   = "true autocorrelation is greater than 0"
  )
  
  names(dw) <- "DW"
  structure(
    list(statistic   = dw,
         method      = "Durbin-Watson test",
         alternative = alternative_str,
         p.value     = pval,
         data.name   = dname),
    class = "htest"
  )
}

.dw_pvalue <- function(dw, X, XtX_inv, n, k,
                       alternative, iterations, exact, tol) {
  if (exact) {
    pval  <- .dw_pvalue_exact(
      dw          = dw,
      X           = X,
      XtX_inv     = XtX_inv,
      n           = n,
      k           = k,
      alternative = alternative,
      iterations  = iterations,
      tol         = tol
    )
    if (is.na(pval) || pval > 1 || pval < 0) {
      warning("exact p value cannot be computed (not in [0,1]), ",
              "approximate p value will be used")
      exact <- FALSE
    }
  }
  
  if (!exact) {
    pval <- .dw_pvalue_approx(
      dw          = dw,
      X           = X,
      XtX_inv     = XtX_inv,
      n           = n,
      k           = k,
      alternative = alternative
    )
  }
  
  pval
}

.dw_pvalue_exact <- function(dw, X, XtX_inv, n, k,
                             alternative, iterations, tol) {
  A        <- diag(c(1, rep(2, n - 2), 1))
  A[abs(row(A) - col(A)) == 1] <- -1
  MA       <- (diag(n) - X %*% XtX_inv %*% t(X)) %*% A
  
  ev_all   <- eigen(MA, only.values = TRUE)$values     # (8) only.values, kein Index-Hack
  
  if (any(abs(Im(ev_all)) > tol))                      # (5) Warnung wieder drin
    warning("imaginary parts of eigenvalues discarded")
  
  ev <- Re(ev_all)
  ev <- ev[ev > tol]
  
  pdw <- function(x) pan(c(x, ev), length(ev), 0, iterations)
  
  switch(alternative,
         "two.sided" = 2 * min(pdw(dw), 1 - pdw(dw)),
         "less"      = 1 - pdw(dw),
         "greater"   = pdw(dw)
  )
}

.dw_pvalue_approx <- function(dw, X, XtX_inv, n, k, alternative) {
  if (n < max(5, k)) {
    warning("not enough observations for computing an approximate p value, set to 1")
    return(1)
  }
  
  AX      <- matrix(as.vector(filter(X, c(-1, 2, -1))), ncol = k)
  AX[1, ] <- X[1, ] - X[2, ]
  AX[n, ] <- X[n, ] - X[n - 1, ]
  
  XAXQ  <- t(X) %*% AX %*% XtX_inv
  P     <- 2 * (n - 1) - sum(diag(XAXQ))
  Q     <- 2 * (3 * n - 4) - 2 * sum(diag(crossprod(AX) %*% XtX_inv)) +
    sum(diag(XAXQ %*% XAXQ))
  dmean <- P / (n - k)
  dvar  <- 2 / ((n - k) * (n - k + 2)) * (Q - P * dmean)
  
  switch(alternative,
         "two.sided" = 2 * pnorm(abs(dw - dmean), sd = sqrt(dvar), lower.tail = FALSE),
         "less"      = pnorm(dw, mean = dmean, sd = sqrt(dvar), lower.tail = FALSE),
         "greater"   = pnorm(dw, mean = dmean, sd = sqrt(dvar))
  )
}


