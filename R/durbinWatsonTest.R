
#' Durbin-Watson Test
#'
#' Tests for first-order autocorrelation in the residuals of a linear
#' regression model, based on the ratio of successive squared residual
#' differences to the total residual sum of squares.
#'
#' The Durbin-Watson test has the null hypothesis that the autocorrelation
#' of the disturbances is 0. The alternative hypothesis can be specified as
#' greater than, not equal to, or less than 0 via the \code{alternative}
#' argument.
#'
#' Under the assumption of normally distributed disturbances, the null
#' distribution of the Durbin-Watson statistic is the distribution of a
#' linear combination of chi-squared variables. The p-value is computed
#' using the "pan" algorithm (Farebrother, 1980, 1984), implemented via
#' Rcpp. For large sample sizes the algorithm might fail to compute the
#' p-value; in that case a warning is issued and an approximate p-value is
#' returned instead, computed via a normal approximation using the mean and
#' variance of the Durbin-Watson statistic.
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
#' @name durbinWatsonTest
#' @param x a symbolic description of the model to be tested (a
#' \code{formula}), a fitted \code{"lm"} object, or a numeric vector of
#' residuals.
#' @param data an optional data frame containing the variables in the
#' model. Only used for the \code{formula} method. By default the variables
#' are taken from the environment which \code{durbinWatsonTest} is called
#' from.
#' @param orderBy either a vector \code{z} or a formula with a single
#' explanatory variable like \code{~ z}. The observations in the model are
#' ordered by the size of \code{z}. If set to \code{NULL} (the default) the
#' observations are assumed to be ordered (e.g., a time series).
#' @param alternative a character string specifying the alternative
#' hypothesis, must be one of \code{"greater"} (default),
#' \code{"two.sided"} or \code{"less"}.
#' @param iterations an integer specifying the number of iterations used by
#' the "pan" algorithm when computing the exact p-value.
#' @param exact logical. If \code{TRUE} the exact p-value is computed via
#' the "pan" algorithm; if \code{FALSE} a normal approximation is used. The
#' default is \code{TRUE} for sample sizes below 100 and \code{FALSE}
#' otherwise.
#' @param tol numeric tolerance. Eigenvalues smaller than \code{tol} are
#' treated as zero.
#' @param \dots further arguments passed to or from other methods.
#'
#' @return An object of class \code{"htest"} containing the following
#' components:
#' \item{statistic}{the Durbin-Watson test statistic.}
#' \item{p.value}{the p-value of the test.}
#' \item{alternative}{a character string describing the alternative
#' hypothesis.}
#' \item{method}{a character string indicating the test performed.}
#' \item{data.name}{a character string describing the data.}
#'
#' @note
#' Based on code by Torsten Hothorn, Achim Zeileis, Clint Cummins,
#' Giovanni Millo and David Mitchell previously published in the
#' \pkg{lmtest} package, with an Rcpp reimplementation of the "pan"
#' algorithm, adapted to conform to package standards.
#'
#' @references
#' Durbin, J. and Watson, G. S. (1950) Testing for serial correlation in
#' least squares regression I. \emph{Biometrika}, 37, 409-428.
#'
#' Durbin, J. and Watson, G. S. (1951) Testing for serial correlation in
#' least squares regression II. \emph{Biometrika}, 38, 159-178.
#'
#' Durbin, J. and Watson, G. S. (1971) Testing for serial correlation in
#' least squares regression III. \emph{Biometrika}, 58, 1-19.
#'
#' Farebrother, R. W. (1980) Pan's procedure for the tail probabilities of
#' the Durbin-Watson statistic. \emph{Applied Statistics}, 29, 224-227.
#'
#' Farebrother, R. W. (1984) AS R53: A remark on algorithms AS 106, AS 153
#' and AS 155. \emph{Applied Statistics}, 33, 366-369.
#'
#' Kraemer, W. and Sonnberger, H. (1986) \emph{The Linear Regression Model
#' under Test}. Heidelberg: Physica.
#'
#' @seealso \code{\link{breuschGodfreyTest}}, \code{\link{lm}}
#'
#' @examples
#' ## formula method
#' set.seed(1)
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
#' ## numeric method (pre-computed residuals, intercept-only design assumed)
#' e_t <- c(-32.33, -26.603, 2.215, -16.967, -1.148, -2.512, -1.967, 11.669,
#'          -0.513, 27.032, -4.422, 40.032, 23.577, 33.94, -2.787, -8.606,
#'           0.575, 6.848, -18.971, -29.063)
#' durbinWatsonTest(e_t)
#'
#' @family test.regression
#' @concept regression-diagnostics
#' @concept autocorrelation
#'
#' @export
durbinWatsonTest <- function(x, ...) UseMethod("durbinWatsonTest")


#' @rdname durbinWatsonTest
#' @export
durbinWatsonTest.formula <- function(x, data = list(), orderBy = NULL,
                                     alternative = c("greater", "two.sided",
                                                     "less"),
                                     iterations = 15, exact = NULL,
                                     tol = 1e-10, ...) {
  mf <- model.frame(x, data = data)
  y  <- model.response(mf)
  X  <- model.matrix(attr(mf, "terms"), mf)

  .dwOrder(X = X, y = y, orderBy = orderBy, data = data,
           dname = deparse1(substitute(x)),
           alternative = alternative, iterations = iterations,
           exact = exact, tol = tol)
}


#' @rdname durbinWatsonTest
#' @export
durbinWatsonTest.lm <- function(x, orderBy = NULL,
                                alternative = c("greater", "two.sided",
                                                "less"),
                                iterations = 15, exact = NULL,
                                tol = 1e-10, ...) {
  if (!is.null(w <- weights(x)))
    if (!isTRUE(all.equal(as.vector(w), rep(1L, length(w)))))
      stop("weighted regressions are not supported")

  X <- if (is.matrix(x$x)) x$x else model.matrix(terms(x), model.frame(x))
  y <- if (is.vector(x$y)) x$y else model.response(model.frame(x))

  .dwOrder(X = X, y = y, orderBy = orderBy, data = list(),
           dname = deparse1(formula(x)),   # the formula, not the object name
           alternative = alternative, iterations = iterations,
           exact = exact, tol = tol)
}


#' @rdname durbinWatsonTest
#' @export
durbinWatsonTest.numeric <- function(x, orderBy = NULL,
                                     alternative = c("greater", "two.sided",
                                                     "less"),
                                     iterations = 15, exact = NULL,
                                     tol = 1e-10, ...) {
  X <- matrix(1, nrow = length(x), ncol = 1)

  .dwOrder(X = X, y = x, orderBy = orderBy, data = list(),
           dname = deparse1(substitute(x)),
           alternative = alternative, iterations = iterations,
           exact = exact, tol = tol)
}


#' @rdname durbinWatsonTest
#' @export
durbinWatsonTest.default <- function(x, ...) {
  stop("no applicable method for objects of class ", sQuote(class(x)[1]))
}



# == internal helper functions ============================================


.dwOrder <- function(X, y, orderBy, data, dname,
                     alternative, iterations, exact, tol) {

  if (!is.null(orderBy)) {
    # accept both a formula ~ z and a plain vector z
    if (inherits(orderBy, "formula")) {
      mm <- model.matrix(orderBy, data = data)
      z  <- as.vector(mm[, ncol(mm)])
    } else {
      z <- orderBy
    }
    X <- as.matrix(X[order(z), ])
    y <- y[order(z)]
  }

  .dwCompute(X = X, y = y, dname = dname, alternative = alternative,
             iterations = iterations, exact = exact, tol = tol)
}


.dwCompute <- function(X, y, dname,
                       alternative = c("greater", "two.sided", "less"),
                       iterations = 15, exact = NULL, tol = 1e-10) {

  alternative <- match.arg(alternative)
  n <- nrow(X)
  k <- ncol(X)
  if (is.null(exact)) exact <- (n < 100)

  fit     <- lm.fit(X, y)
  res     <- fit$residuals
  dw      <- sum(diff(res)^2) / sum(res^2)
  XtX_inv <- chol2inv(qr.R(fit$qr))

  pval <- if (n < 3) {
    warning("not enough observations for computing a p value, set to 1")
    1
  } else {
    .dwPvalue(dw = dw, X = X, XtX_inv = XtX_inv, n = n, k = k,
              alternative = alternative, iterations = iterations,
              exact = exact, tol = tol)
  }

  ALTERNATIVE <- switch(alternative,
                        "two.sided" = "true autocorrelation is not 0",
                        "less"      = "true autocorrelation is less than 0",
                        "greater"   = "true autocorrelation is greater than 0")

  structure(
    list(statistic   = c(DW = dw),
         method      = "Durbin-Watson test",
         alternative = ALTERNATIVE,
         p.value     = pval,
         data.name   = dname),
    class = "htest")
}


.dwPvalue <- function(dw, X, XtX_inv, n, k,
                      alternative, iterations, exact, tol) {

  if (exact) {
    pval <- .dwPvalueExact(dw = dw, X = X, XtX_inv = XtX_inv, n = n, k = k,
                           alternative = alternative,
                           iterations = iterations, tol = tol)
    if (is.na(pval) || pval > 1 || pval < 0) {
      warning("exact p value cannot be computed (not in [0,1]), ",
              "approximate p value will be used")
      exact <- FALSE
    }
  }

  if (!exact)
    pval <- .dwPvalueApprox(dw = dw, X = X, XtX_inv = XtX_inv, n = n, k = k,
                            alternative = alternative)

  pval
}


.dwPvalueExact <- function(dw, X, XtX_inv, n, k,
                           alternative, iterations, tol) {

  A <- diag(c(1, rep(2, n - 2), 1))
  A[abs(row(A) - col(A)) == 1] <- -1
  MA <- (diag(n) - X %*% XtX_inv %*% t(X)) %*% A

  ev_all <- eigen(MA, only.values = TRUE)$values

  if (any(abs(Im(ev_all)) > tol))
    warning("imaginary parts of eigenvalues discarded")

  ev <- Re(ev_all)
  ev <- ev[ev > tol]

  # pan() returns P(DW <= x), see the Rcpp implementation
  pdw <- function(x) pan_cpp(c(x, ev), length(ev), 0, iterations)

  switch(alternative,
         "two.sided" = 2 * min(pdw(dw), 1 - pdw(dw)),
         "less"      = 1 - pdw(dw),
         "greater"   = pdw(dw))
}



.dwPvalueApprox <- function(dw, X, XtX_inv, n, k, alternative) {

  if (n < max(5, k)) {
    warning("not enough observations for computing an approximate p value, ",
            "set to 1")
    return(1)
  }

  AX      <- matrix(as.vector(stats::filter(X, c(-1, 2, -1))), ncol = k)
  AX[1, ] <- X[1, ] - X[2, ]
  AX[n, ] <- X[n, ] - X[n - 1, ]

  XAXQ  <- t(X) %*% AX %*% XtX_inv
  P     <- 2 * (n - 1) - sum(diag(XAXQ))
  Q     <- 2 * (3 * n - 4) - 2 * sum(diag(crossprod(AX) %*% XtX_inv)) +
    sum(diag(XAXQ %*% XAXQ))
  dmean <- P / (n - k)
  dvar  <- 2 / ((n - k) * (n - k + 2)) * (Q - P * dmean)

  switch(alternative,
         "two.sided" = 2 * pnorm(abs(dw - dmean), sd = sqrt(dvar),
                                 lower.tail = FALSE),
         "less"      = pnorm(dw, mean = dmean, sd = sqrt(dvar),
                             lower.tail = FALSE),
         "greater"   = pnorm(dw, mean = dmean, sd = sqrt(dvar)))
}
