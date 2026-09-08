
#' Breusch-Godfrey Test for Detecting Higher-Order Serial Correlation in Regression Residuals
#'
#' A test for autocorrelation in the residuals of regression models,
#' generalizing the Durbin-Watson test to handle higher-order
#' autocorrelation and models with lagged dependent variables.
#'
#' `breuschGodfreyTest` performs the Breusch-Godfrey test for
#' higher-order serial correlation.
#'
#' Under \eqn{H_0} the test statistic is asymptotically chi-squared with
#' degrees of freedom as given in `parameter`. If `type` is set
#' to `"f"` the function returns a finite sample version of the test
#' statistic, employing an \eqn{F} distribution with degrees of freedom as
#' given in `parameter`.
#'
#' By default, the starting values for the lagged residuals in the
#' auxiliary regression are chosen to be 0 (as in Godfrey 1978) but could
#' also be set to `NA` to omit them.
#'
#' `breuschGodfreyTest` also returns the coefficients and estimated
#' covariance matrix from the auxiliary regression that includes the lagged
#' residuals, accessible via `coef()` and `vcov()` on the result.
#' (Note, however, that standard theory does not always apply to the
#' standard errors and t-statistics in this regression.)
#'
#' @param formula a symbolic description for the model to be tested (or a
#' fitted `"lm"` object, in which case the model frame is taken from the fit
#' and `subset` and `na.action` are ignored).
#' @param data an optional data frame containing the variables in the
#' model. By default the variables are taken from the environment which
#' `breuschGodfreyTest` is called from. For a fitted `"lm"` object it is used
#' for `orderBy` only, as the model itself carries its own model frame.
#' @param order integer, the maximal order of serial correlation to be
#' tested. Must be smaller than the residual degrees of freedom of the
#' auxiliary regression.
#' @param orderBy either a vector `z` or a formula with a single
#' explanatory variable like `~ z`. The observations in the model are
#' ordered by the size of `z`; a formula with several terms is used as
#' successive ordering keys. If set to `NULL` (the default) the observations
#' are assumed to be ordered (e.g., a time series). `z` may be given at the
#' length of the original data: rows dropped by `subset` or by `na.action`
#' are then dropped from `z` as well. Missing values in `z` are ordered
#' last.
#' @param type the type of test statistic to be returned, either
#' `"chisq"` (default) for the chi-squared test statistic or
#' `"f"` for the F test statistic. Case-insensitive.
#' @param subset an optional expression indicating which observations to use.
#' @param na.action a function specifying how missing values are handled.
#' Defaults to [na.omit()]: the auxiliary regression is fitted by
#' [lm.fit()] and cannot carry missing values.
#' @param fill a single value used as starting value for the lagged residuals
#' in the auxiliary regression. By default `0` but can also be set to `NA`,
#' in which case the leading incomplete rows are dropped.
#'
#' @return A list with class `"breuschGodfreyTest"` inheriting from
#' `"htest"` containing the following components:
#'   \item{`statistic`}{the value of the test statistic.}
#'   \item{`parameter`}{the degrees of freedom.}
#'   \item{`p.value`}{the p-value of the test.}
#'   \item{`method`}{a character string indicating what type of test was
#'     performed.}
#'   \item{`data.name`}{a character string giving the name(s) of the
#'     data.}
#'   \item{`coefficients`}{coefficient estimates from the auxiliary
#'     regression.}
#'   \item{`vcov`}{the corresponding covariance matrix estimate.}
#'   \item{`df.residual`}{the residual degrees of freedom of the auxiliary
#'     regression, for both types of test statistic.}
#'
#' @note
#' Based on code by David Mitchell and Achim Zeileis previously published
#' as `bgtest()` in the \pkg{lmtest} package, adapted to conform to
#' package standards.
#'
#' Unlike `bgtest()`, the residual degrees of freedom of the auxiliary
#' regression are reported for both types of test statistic. They are a
#' property of that regression and not of the statistic derived from it,
#' whereas `bgtest()` hands back `NULL` for the chi-squared version.
#' `coeftest()` therefore refers the coefficients to a \eqn{t} distribution
#' in either case, where `bgtest()` switches to the normal one.
#'
#' @references
#' Breusch, T. S. (1978) Testing for autocorrelation in dynamic linear
#' models. *Australian Economic Papers*, 17, 334-355.
#'
#' Godfrey, L. G. (1978) Testing against general autoregressive and moving
#' average error models when the regressors include lagged dependent
#' variables. *Econometrica*, 46, 1293-1301.
#'
#' @examples
#' ## Generate a stationary and an AR(1) series
#' set.seed(1)
#' x <- rep(c(1, -1), 50)
#'
#' y1 <- 1 + x + rnorm(100)
#'
#' ## Perform Breusch-Godfrey test for first-order serial correlation:
#' breuschGodfreyTest(y1 ~ x)
#'
#' ## or for fourth-order serial correlation
#' breuschGodfreyTest(y1 ~ x, order = 4)
#'
#' ## Compare with Durbin-Watson test results:
#' durbinWatsonTest(y1 ~ x)
#'
#' y2 <- stats::filter(y1, 0.5, method = "recursive")
#' breuschGodfreyTest(y2 ~ x)
#'
#' ## finite sample F version, and dropping the leading lags instead of
#' ## filling them with zeros
#' breuschGodfreyTest(y2 ~ x, order = 4, type = "f")
#' breuschGodfreyTest(y2 ~ x, order = 4, fill = NA)
#'
#' ## transformed terms and an explicit ordering variable
#' d <- data.frame(y = as.vector(y2), x = x, z = rnorm(100), tt = sample(100),
#'                 grp = rep(c("A", "B"), each = 50))
#' breuschGodfreyTest(y ~ x + I(z^2), data = d, orderBy = ~ tt)
#'
#' ## subset and orderBy combined: tt is given at the length of d and is
#' ## reduced to the rows the model frame kept
#' breuschGodfreyTest(y ~ x, data = d, subset = grp == "A", orderBy = ~ tt)
#'
#' ## the test can also be applied to a fitted model
#' breuschGodfreyTest(lm(y1 ~ x))
#'
#' @seealso [durbinWatsonTest()]
#'
#' @family test.regression
#' @concept regression-diagnostics
#' @concept autocorrelation
#'
#' @export
breuschGodfreyTest <- function(formula, data = list(), order = 1,
                               orderBy = NULL, type = c("chisq", "f"),
                               subset, na.action = na.omit, fill = 0) {

  # ── Validate ──────────────────────────────────────────────────────────────
  type <- match.arg(tolower(type), c("chisq", "f"))

  nLags <- suppressWarnings(as.integer(order))
  if (length(nLags) != 1L || is.na(nLags) || nLags < 1L)
    stop("'order' must be a single positive integer", call. = FALSE)

  if (length(fill) != 1L ||
      !(is.numeric(fill) || (is.logical(fill) && is.na(fill))))
    stop("'fill' must be a single numeric value or NA", call. = FALSE)

  # ── Response and design matrix ────────────────────────────────────────────
  if (inherits(formula, "formula")) {

    subsetExpr <- if (missing(subset)) NULL else substitute(subset)

    r <- resolveFormula(formula, data = data, subset = subsetExpr,
                        na.action = na.action, allowed = "regression")

    y <- r$response
    # the model matrix must be built from the terms: the columns of a model
    # frame are named after the deparsed expressions ("log(x)"), so the
    # formula itself cannot be re-evaluated against it
    X     <- model.matrix(r$terms, r$mf)
    rows  <- r$rows
    dname <- r$dataName

  } else {

    dname <- deparse1(substitute(formula))
    # [[exact = TRUE]] rather than $: partial matching would resolve $x to
    # the xlevels component of an "lm" object
    xComp <- formula[["x", exact = TRUE]]
    yComp <- formula[["y", exact = TRUE]]

    X <- if (is.matrix(xComp)) xComp
         else model.matrix(terms(formula), model.frame(formula))
    y <- if (is.vector(yComp)) yComp
         else model.response(model.frame(formula))

    # the counterpart of resolveFormula()'s 'rows' for a fitted model
    rows <- .rowsFromNames(rownames(X), data)
  }

  if (NCOL(y) != 1L)
    stop("the response must be a single vector", call. = FALSE)

  y <- as.vector(y)

  if (is.null(colnames(X)))
    colnames(X) <- paste0("x", seq_len(ncol(X)))

  # ── Reorder ───────────────────────────────────────────────────────────────
  ord <- .orderIndex(orderBy, nrow(X), data = data, rows = rows)

  if (!is.null(ord)) {
    X <- X[ord, , drop = FALSE]
    y <- y[ord]
  }

  # ── Auxiliary regression ──────────────────────────────────────────────────
  n <- nrow(X)
  k <- ncol(X)

  if (nLags >= n)
    stop("'order' must be smaller than the number of observations",
         call. = FALSE)

  fit <- lm.fit(X, y)

  # a defect of the model itself, reported here rather than as a rank
  # deficiency of the auxiliary regression below
  if (fit$rank < k)
    stop("the model matrix is rank deficient", call. = FALSE)

  resi <- fit$residuals
  lags <- seq_len(nLags)

  Z <- vapply(lags, function(i) c(rep(fill, i), resi[seq_len(n - i)]),
              numeric(n))

  if (any(incomplete <- !complete.cases(Z))) {
    X    <- X[!incomplete, , drop = FALSE]
    Z    <- Z[!incomplete, , drop = FALSE]
    resi <- resi[!incomplete]
    n    <- nrow(X)
  }

  if (n - k - nLags < 1L)
    stop(gettextf("not enough observations to test up to order %d", nLags),
         call. = FALSE)

  auxfit <- lm.fit(cbind(X, Z), resi)

  if (auxfit$rank < k + nLags)
    stop("the auxiliary regression is rank deficient", call. = FALSE)

  cf <- auxfit$coefficients
  vc <- chol2inv(auxfit$qr$qr) *
    sum(auxfit$residuals^2) / auxfit$df.residual
  names(cf) <- colnames(vc) <- rownames(vc) <-
    c(colnames(X), paste("lag(resid)", lags, sep = "_"))

  # ── Statistic ─────────────────────────────────────────────────────────────
  switch(type,

         "chisq" = {
           bg    <- n * sum(auxfit$fitted.values^2) / sum(resi^2)
           df    <- c(df = nLags)
           p.val <- unname(pchisq(bg, nLags, lower.tail = FALSE))
         },

         "f" = {
           ssrU  <- sum(auxfit$residuals^2)
           bg    <- ((sum(resi^2) - ssrU) / nLags) / (ssrU / auxfit$df.residual)
           df    <- c(df1 = nLags, df2 = auxfit$df.residual)
           p.val <- unname(pf(bg, df1 = df[1L], df2 = df[2L],
                              lower.tail = FALSE))
         })

  names(bg) <- "LM test"

  method <- gettextf(
    "Breusch-Godfrey test for serial correlation of order up to %d", nLags)

  res <- list(statistic   = bg,
              parameter   = df,
              method      = method,
              p.value     = p.val,
              data.name   = dname,
              coefficients = cf,
              vcov        = vc,
              df.residual = auxfit$df.residual)

  class(res) <- c("breuschGodfreyTest", "htest")
  res
}


#' @export
vcov.breuschGodfreyTest <- function(object, ...)
  object$vcov


#' @export
df.residual.breuschGodfreyTest <- function(object, ...)
  object$df.residual
