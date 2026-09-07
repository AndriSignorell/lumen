
#' Durbin-Watson Test for Detecting First-Order Autocorrelation in Regression Residuals
#'
#' Tests for first-order autocorrelation in the residuals of a linear
#' regression model, based on the ratio of successive squared residual
#' differences to the total residual sum of squares.
#'
#' The Durbin-Watson test has the null hypothesis that the autocorrelation
#' of the disturbances is 0. The alternative hypothesis can be specified as
#' greater than, not equal to, or less than 0 via the `alternative`
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
#'   \item{`formula`}{Fits the model from scratch using `data`.}
#'   \item{`lm`}{Extracts design matrix and response from a fitted
#'     `"lm"` object. The `data.name` field in the result reflects
#'     the model formula, not the name of the object.}
#'   \item{`numeric`}{Treats `x` as a vector of pre-computed
#'     residuals and fits an intercept-only design matrix. Note that this is
#'     *not* equivalent to testing residuals from a fitted model with
#'     predictors; use the `lm` or `formula` method when a model
#'     exists.}
#' }
#'
#' @param x a symbolic description of the model to be tested (a
#' `formula`), a fitted `"lm"` object, or a numeric vector of
#' residuals.
#' @param data an optional data frame containing the variables in the
#' model. By default the variables are taken from the environment which
#' `durbinWatsonTest` is called from. For the `lm` and `numeric` methods it
#' is used for `orderBy` only, as the model frame is already fixed there.
#' @param orderBy either a vector `z` or a one-sided formula like `~ z`. The
#' observations in the model are ordered by the size of `z`; a formula with
#' several terms is used as successive ordering keys. If set to `NULL` (the
#' default) the observations are assumed to be ordered (e.g., a time
#' series). `z` may be given at the length of the original data: rows
#' dropped by `subset` or by `na.action` are then dropped from `z` as well.
#' Missing values in `z` are ordered last.
#' @param subset an optional expression indicating which observations to
#' use. Only used for the `formula` method.
#' @param na.action a function specifying how missing values are handled.
#' Defaults to [na.omit()]. Only used for the `formula` method.
#' @param alternative a character string specifying the alternative
#' hypothesis, must be one of `"greater"` (default),
#' `"two.sided"` or `"less"`.
#' @param iterations an integer specifying the number of iterations used by
#' the "pan" algorithm when computing the exact p-value.
#' @param exact logical. If `TRUE` the exact p-value is computed via
#' the "pan" algorithm; if `FALSE` a normal approximation is used. The
#' default is `TRUE` for sample sizes below 100 and `FALSE`
#' otherwise.
#' @param tol numeric tolerance. Eigenvalues smaller than `tol` are
#' treated as zero.
#' @param \dots further arguments passed to or from other methods.
#'
#' @return An object of class `"htest"` containing the following
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
#' least squares regression I. *Biometrika*, 37, 409-428.
#'
#' Durbin, J. and Watson, G. S. (1951) Testing for serial correlation in
#' least squares regression II. *Biometrika*, 38, 159-178.
#'
#' Durbin, J. and Watson, G. S. (1971) Testing for serial correlation in
#' least squares regression III. *Biometrika*, 58, 1-19.
#'
#' Farebrother, R. W. (1980) Pan's procedure for the tail probabilities of
#' the Durbin-Watson statistic. *Applied Statistics*, 29, 224-227.
#'
#' Farebrother, R. W. (1984) AS R53: A remark on algorithms AS 106, AS 153
#' and AS 155. *Applied Statistics*, 33, 366-369.
#'
#' Kraemer, W. and Sonnberger, H. (1986) *The Linear Regression Model
#' under Test*. Heidelberg: Physica.
#'
#' @seealso [lm()], [breuschGodfreyTest()]
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
#' ## subset and an ordering variable given at the length of the data
#' d <- data.frame(y = 1 + x + as.vector(err2), x = x, tt = sample(100),
#'                 grp = rep(c("A", "B"), each = 50))
#' durbinWatsonTest(y ~ x, data = d, subset = grp == "A", orderBy = ~ tt)
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
                                     tol = 1e-10,
                                     subset, na.action = na.omit, ...) {

  subsetExpr <- if (missing(subset)) NULL else substitute(subset)

  r <- resolveFormula(x, data = data, subset = subsetExpr,
                      na.action = na.action, allowed = "regression")

  # the model matrix must be built from the terms: the columns of a model
  # frame are named after the deparsed expressions ("log(x)"), so the
  # formula itself cannot be re-evaluated against it
  .dwCompute(X = model.matrix(r$terms, r$mf), y = r$response,
             orderBy = orderBy, data = data, rows = r$rows,
             dname = r$dataName,
             alternative = alternative, iterations = iterations,
             exact = exact, tol = tol)
}


#' @rdname durbinWatsonTest
#' @export
durbinWatsonTest.lm <- function(x, data = list(), orderBy = NULL,
                                alternative = c("greater", "two.sided",
                                                "less"),
                                iterations = 15, exact = NULL,
                                tol = 1e-10, ...) {

  if (!is.null(w <- weights(x)))
    if (!isTRUE(all.equal(as.vector(w), rep(1, length(w)))))
      stop("weighted regressions are not supported", call. = FALSE)

  # [[exact = TRUE]] rather than $: partial matching would resolve $x to the
  # xlevels component of an "lm" object
  xComp <- x[["x", exact = TRUE]]
  yComp <- x[["y", exact = TRUE]]

  X <- if (is.matrix(xComp)) xComp
       else model.matrix(terms(x), model.frame(x))
  y <- if (is.vector(yComp)) yComp
       else model.response(model.frame(x))

  .dwCompute(X = X, y = y,
             orderBy = orderBy, data = data,
             rows = .rowsFromNames(rownames(X), data),
             dname = deparse1(formula(x)),   # the formula, not the object name
             alternative = alternative, iterations = iterations,
             exact = exact, tol = tol)
}


#' @rdname durbinWatsonTest
#' @export
durbinWatsonTest.numeric <- function(x, data = list(), orderBy = NULL,
                                     alternative = c("greater", "two.sided",
                                                     "less"),
                                     iterations = 15, exact = NULL,
                                     tol = 1e-10, ...) {

  # the implicit class of a double matrix contains "numeric", so this method
  # is also reached with a matrix, where length() would build a design of the
  # wrong size
  if (!is.vector(x) || !is.numeric(x))
    stop("'x' must be a numeric vector of residuals", call. = FALSE)

  .dwCompute(X = matrix(1, nrow = length(x), ncol = 1L), y = x,
             orderBy = orderBy, data = data,
             rows = .rowsFromNames(names(x), data),
             dname = deparse1(substitute(x)),
             alternative = alternative, iterations = iterations,
             exact = exact, tol = tol)
}


#' @rdname durbinWatsonTest
#' @export
durbinWatsonTest.default <- function(x, ...) {
  stop(gettextf("no applicable method for objects of class %s",
                sQuote(class(x)[1L])), call. = FALSE)
}



# == internal helper functions ============================================


.dwCompute <- function(X, y, orderBy, data, rows, dname,
                       alternative = c("greater", "two.sided", "less"),
                       iterations = 15, exact = NULL, tol = 1e-10) {

  # ── Validate ──────────────────────────────────────────────────────────────
  # one place for all three methods, before any work is done
  alternative <- match.arg(alternative)

  iterations <- suppressWarnings(as.integer(iterations))
  if (length(iterations) != 1L || is.na(iterations) || iterations < 1L)
    stop("'iterations' must be a single positive integer", call. = FALSE)

  if (length(tol) != 1L || !is.numeric(tol) || is.na(tol) || tol < 0)
    stop("'tol' must be a single non-negative number", call. = FALSE)

  if (!is.null(exact) &&
      (length(exact) != 1L || !is.logical(exact) || is.na(exact)))
    stop("'exact' must be a single logical value or NULL", call. = FALSE)

  if (NCOL(y) != 1L)
    stop("the response must be a single vector", call. = FALSE)

  y <- as.vector(y)

  # ── Reorder ───────────────────────────────────────────────────────────────
  ord <- .orderIndex(orderBy, nrow(X), data = data, rows = rows)

  if (!is.null(ord)) {
    X <- X[ord, , drop = FALSE]
    y <- y[ord]
  }

  # ── Statistic ─────────────────────────────────────────────────────────────
  n <- nrow(X)
  k <- ncol(X)

  # below that the residuals are identically zero and the statistic is 0/0
  if (n < 3L || n <= k)
    stop(gettextf("at least %d observations are needed to compute the statistic",
                  max(3L, k + 1L)), call. = FALSE)

  if (is.null(exact))
    exact <- (n < 100L)

  fit <- lm.fit(X, y)

  if (fit$rank < k)
    stop("the design matrix is rank deficient", call. = FALSE)

  res     <- fit$residuals
  dw      <- sum(diff(res)^2) / sum(res^2)
  XtX_inv <- chol2inv(qr.R(fit$qr))

  pval <- .dwPvalue(dw = dw, X = X, XtX_inv = XtX_inv, n = n, k = k,
                    alternative = alternative, iterations = iterations,
                    exact = exact, tol = tol)

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
    pval <- .dwPvalueExact(dw = dw, X = X, XtX_inv = XtX_inv, n = n,
                           alternative = alternative,
                           iterations = iterations, tol = tol)

    # a failing pan() may return nothing at all, not just something outside
    # the unit interval
    if (length(pval) != 1L || is.na(pval) || pval > 1 || pval < 0) {
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


.dwPvalueExact <- function(dw, X, XtX_inv, n, alternative, iterations, tol) {

  # A is the second difference matrix, M = I - X (X'X)^-1 X' the residual
  # maker. The non-zero eigenvalues of MA are those of the symmetric MAM,
  # which is what the null distribution of DW is built from, so the
  # eigenvalues come out real by construction and need no cleaning up.
  A   <- diag(c(1, rep(2, n - 2), 1))
  idx <- cbind(seq_len(n - 1L), 2:n)
  A[idx] <- A[idx[, 2:1]] <- -1

  # M %*% A %*% M without ever forming M: O(n^2 k) instead of O(n^3)
  MA  <- A - X %*% (XtX_inv %*% crossprod(X, A))
  MAM <- MA - tcrossprod(MA %*% X %*% XtX_inv, X)

  ev <- eigen(MAM, symmetric = TRUE, only.values = TRUE)$values
  ev <- ev[ev > tol]

  # pan() returns P(DW <= x), see the Rcpp implementation
  p <- pan_cpp(c(dw, ev), length(ev), 0, iterations)

  switch(alternative,
         "two.sided" = 2 * min(p, 1 - p),
         "less"      = 1 - p,
         "greater"   = p)
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

  if (!is.finite(dvar) || dvar <= 0) {
    warning("the variance of the statistic is not positive, ",
            "approximate p value set to 1")
    return(1)
  }

  switch(alternative,
         "two.sided" = 2 * pnorm(abs(dw - dmean), sd = sqrt(dvar),
                                 lower.tail = FALSE),
         "less"      = pnorm(dw, mean = dmean, sd = sqrt(dvar),
                             lower.tail = FALSE),
         "greater"   = pnorm(dw, mean = dmean, sd = sqrt(dvar)))
}
