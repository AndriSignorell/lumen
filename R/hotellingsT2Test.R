#' Hotelling's T2 Test
#'
#' The classical test for the location of a multivariate population (one-sample)
#' or for the difference between the mean vectors of two multivariate populations
#' (two-sample).  It is the multivariate generalisation of Student's t-test.
#'
#' @details
#' When \code{test = "f"} the test statistic follows an exact F-distribution
#' under the assumption of multivariate normality.  When \code{test = "chi"} a
#' chi-squared approximation is used; it relies on large-sample asymptotic
#' theory and is less sensitive to departures from multivariate normality than
#' the F-test, but remains only asymptotically correct.
#'
#' In the two-sample case both populations are assumed to share the same
#' covariance matrix; a pooled within-group estimate is used.
#'
#' The formula interface (\code{cbind(v1, v2) ~ g}) is available for the
#' two-sample case only.
#'
#' @param x       A numeric matrix or data frame (observations in rows,
#'   variables in columns).
#' @param y       An optional numeric matrix or data frame for the two-sample
#'   test.  If \code{NULL} (default) a one-sample test is performed.
#' @param mu      A numeric vector of length \eqn{p} giving the hypothesised
#'   mean (one-sample) or mean difference (two-sample).  \code{NULL} is
#'   interpreted as the zero vector.
#' @param test    A character string selecting the reference distribution:
#'   \code{"f"} (exact F-distribution, default) or \code{"chi"}
#'   (chi-squared approximation).
#' @param formula A formula of the form \code{cbind(v1, v2, ...) ~ g} where
#'   the left-hand side is a numeric matrix of response variables and \code{g}
#'   is a factor with exactly two levels.
#' @param data    An optional data frame (or similar, see
#'   \code{\link{model.frame}}) containing the variables in \code{formula}.
#'   Defaults to the environment of \code{formula}.
#' @param subset  An optional vector specifying a subset of observations.
#' @param na.action A function indicating what should happen when the data
#'   contain \code{NA}s.  Defaults to \code{getOption("na.action")}.
#' @param \dots  Further arguments passed to or from methods.
#'
#' @return An object of class \code{"htest"} containing:
#'   \item{statistic}{The value of the T2-statistic (scaled to follow an F- or
#'     chi-squared distribution depending on \code{test}).}
#'   \item{parameter}{Degrees of freedom of the reference distribution.}
#'   \item{p.value}{The p-value of the test.}
#'   \item{null.value}{The hypothesised mean or mean difference.}
#'   \item{alternative}{Always \code{"two.sided"}.}
#'   \item{method}{A character string describing the test variant performed.}
#'   \item{data.name}{A character string giving the name(s) of the input
#'     data.}
#'
#' @note
#' Based on code by Klaus Nordhausen, adapted to conform to package standards.
#'
#' @references
#' Anderson, T. W. (2003). \emph{An Introduction to Multivariate Statistical
#'   Analysis} (3rd ed.). Wiley.
#'
#' Nordhausen, K., Sirkia, S., Oja, H., & Tyler, D. E. (2012).
#'   \emph{ICSNP: Tools for Multivariate Nonparametrics}. R package
#'   version 1.0-9. \url{https://cran.r-project.org/package=ICSNP}
#'
#' @examples
#' math.teach <- data.frame(
#'   teacher = factor(rep(1:2, c(3, 6))),
#'   satis   = c(1, 3, 2, 4, 6, 6, 5, 5, 4),
#'   know    = c(3, 7, 2, 6, 8, 8, 10, 10, 6))
#'
#' with(math.teach,
#'   hotellingsT2Test(cbind(satis, know) ~ teacher))
#'
#' # chi-squared approximation
#' with(math.teach,
#'   hotellingsT2Test(cbind(satis, know) ~ teacher, test = "chi"))
#'
#' @name hotelling-t2-test
#' @aliases hotellingsT2Test hotellingsT2Test.default hotellingsT2Test.formula
#' @family topic.hypothesisTests
#' @concept hypothesis-testing
#' @concept multivariate
#' @concept location
#'
#' @export
hotellingsT2Test <- function(x, ...) {
  UseMethod("hotellingsT2Test")
}



#' @rdname hotelling-t2-test
#' @export
hotellingsT2Test.formula <- function(formula,
                                     data,
                                     subset,
                                     na.action = na.pass,
                                     ...) {
  
  if (missing(formula) || length(formula) != 3L)
    stop("'formula' missing or incorrect.")
  
  args <- list(
    formula   = formula,
    na.action = na.action,
    allowed   = "two.sample.independent"
  )
  if (!missing(data))   args$data   <- data
  if (!missing(subset)) args$subset <- substitute(subset)
  
  rf <- do.call(bedrock::resolveFormula, args)
  
  # resolveFormula uses split() internally, which flattens a matrix response
  # to a vector.  Restore the original column structure.
  p <- ncol(rf$mf[[1L]])
  x <- matrix(rf$x, ncol = p, dimnames = list(NULL, colnames(rf$mf[[1L]])))
  y <- matrix(rf$y, ncol = p, dimnames = list(NULL, colnames(rf$mf[[1L]])))
  
  res           <- hotellingsT2Test.default(x = x, y = y, ...)
  res$data.name <- rf$data.name
  res
}



#' @rdname hotelling-t2-test
#' @export
hotellingsT2Test.default <- function(x, y = NULL, mu = NULL,
                                     test = c("f", "chi"), ...) {
  
  # --- data name -------------------------------------------------------
  DNAME <- if (is.null(y)) {
    deparse(substitute(x))
  } else {
    paste(deparse(substitute(x)), "and", deparse(substitute(y)))
  }
  
  # --- validate x ------------------------------------------------------
  x <- as.matrix(x)
  if (!is.numeric(x))
    stop("Argument 'x' must be numeric.")
  x <- x[complete.cases(x), , drop = FALSE]
  p <- ncol(x)
  
  if (nrow(x) <= p)
    stop("'x' must have more rows than columns for the covariance matrix ",
         "to be non-singular.")
  
  # --- validate y ------------------------------------------------------
  if (!is.null(y)) {
    y <- as.matrix(y)
    if (!is.numeric(y))
      stop("Argument 'y' must be numeric.")
    y <- y[complete.cases(y), , drop = FALSE]
    if (ncol(y) != p)
      stop("'x' and 'y' must have the same number of columns.")
    
    if (nrow(y) <= p)
      stop("'y' must have more rows than columns for the covariance matrix ",
           "to be non-singular.")
  }
  
  # --- validate mu -----------------------------------------------------
  if (is.null(mu)) {
    mu <- rep(0, p)
  } else {
    if (!is.numeric(mu) || !all(is.finite(mu)))
      stop("Argument 'mu' must be a finite numeric vector.")
    if (length(mu) != p)
      stop("Length of 'mu' must equal the number of columns of 'x' (", p, ").")
  }
  
  test <- match.arg(test)
  
  # --- compute ---------------------------------------------------------
  res <- .hotellingsT2_engine(x, y, mu, test)
  
  # --- build htest object ----------------------------------------------
  STATISTIC           <- res$statistic
  names(STATISTIC)    <- "T.2"
  
  PARAMETER <- if (!is.na(res$df2)) {
    c(df1 = res$df1, df2 = res$df2)
  } else {
    c(df = res$df1)
  }
  
  METHOD <- if (is.null(y)) {
    "Hotelling's one-sample T-squared test"
  } else {
    "Hotelling's two-sample T-squared test"
  }
  
  nval_label  <- if (is.null(y)) "location" else "location difference"
  NVAL        <- setNames(mu, rep(nval_label, length(mu)))
  
  structure(
    list(
      statistic   = STATISTIC,
      parameter   = PARAMETER,
      p.value     = res$pValue,
      null.value  = NVAL,
      alternative = "two.sided",
      method      = METHOD,
      data.name   = DNAME
    ),
    class = "htest"
  )
}



# == internal helper functions =================================================


# Internal engine — operates on validated matrices only

.hotellingsT2_engine <- function(x, y = NULL, mu, test) {
  
  n <- nrow(x)
  p <- ncol(x)
  
  if (is.null(y)) {
    
    # one-sample
    diff    <- colMeans(x) - mu
    S_inv   <- tryCatch(
      chol2inv(chol(cov(x))),
      error = function(e) stop(
        "Covariance matrix of 'x' is singular or nearly singular; ",
        "the test cannot be computed."
      )
    )
    T2      <- as.numeric(n * t(diff) %*% S_inv %*% diff)
    scale   <- switch(test, f = (n - p) / (p * (n - 1)), chi = 1)
    stat    <- T2 * scale
    df1     <- p
    df2     <- switch(test, f = n - p,      chi = NA_real_)
    pval    <- switch(test,
                      f   = pf(stat,        df1, df2, lower.tail = FALSE),
                      chi = pchisq(stat,    df1,       lower.tail = FALSE))
    
  } else {
    
    # two-sample
    n1      <- n
    n2      <- nrow(y)
    xm      <- colMeans(x)
    ym      <- colMeans(y)
    Sp      <- (crossprod(sweep(x, 2, xm)) +
                  crossprod(sweep(y, 2, ym))) / (n1 + n2 - 2)
    Sp_inv  <- tryCatch(
      chol2inv(chol(Sp)),
      error = function(e) stop(
        "Pooled covariance matrix is singular or nearly singular; ",
        "the test cannot be computed."
      )
    )
    diff    <- xm - ym - mu
    T2      <- as.numeric(n1 * n2 / (n1 + n2) * t(diff) %*% Sp_inv %*% diff)
    scale   <- switch(test,
                      f   = (n1 + n2 - p - 1) / (p * (n1 + n2 - 2)),
                      chi = 1)
    stat    <- T2 * scale
    df1     <- p
    df2     <- switch(test, f = n1 + n2 - p - 1, chi = NA_real_)
    pval    <- switch(test,
                      f   = pf(stat,        df1, df2, lower.tail = FALSE),
                      chi = pchisq(stat,    df1,       lower.tail = FALSE))
  }
  
  list(statistic = stat, pValue = pval, df1 = df1, df2 = df2)
}

