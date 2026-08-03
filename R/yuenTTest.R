
#' Yuen T-Test for Robust Comparison of Trimmed Means
#'
#' Robust one-, two-, and paired-sample t-tests based on trimmed means
#' and winsorized variances.
#'
#' @description
#' Performs Yuen's robust t-test for trimmed means. Compared with the
#' classical t-test, the procedure is substantially less sensitive to
#' outliers, heavy tails, and moderate departures from normality.
#'
#' The test is based on:
#' \itemize{
#'   \item trimmed means,
#'   \item winsorized variances,
#'   \item Welch-type degrees of freedom.
#' }
#'
#' For paired tests, trimming is performed on the paired differences,
#' following Yuen (1974).
#'
#' @name yuenTTest
#' @aliases yuenTTest yuenTTest.default yuenTTest.formula
#'
#' @param x numeric vector of observations.
#' @param y optional second numeric vector.
#' @param alternative character string specifying the alternative
#'   hypothesis. One of \code{"two.sided"}, \code{"less"},
#'   or \code{"greater"}.
#' @param paired logical indicating whether a paired test is performed.
#' @param mu hypothesized trimmed mean (or trimmed mean difference).
#' @param conf.level confidence level for the confidence interval.
#' @param trim fraction of observations trimmed from each tail.
#'   Must satisfy \code{0 <= trim < 0.5}.
#' @param formula a formula of the form \code{lhs ~ rhs}.
#' @param data optional data frame for the formula interface.
#' @param subset optional subset expression.
#' @param na.action NA handling function.
#' @param \dots further arguments passed to methods.
#'
#' @return
#' An object of class \code{"htest"}.
#'
#' @seealso
#' \code{\link{t.test}}
#'
#' @references
#' Wilcox, R. R. (2005).
#' \emph{Introduction to Robust Estimation and Hypothesis Testing}.
#' Academic Press.
#'
#' Yuen, K. K. (1974).
#' The two-sample trimmed t for unequal population variances.
#' \emph{Biometrika}, 61, 165--170.
#'
#' @examples
#' x <- rnorm(25, 100, 5)
#' yuenTTest(x, mu = 99)
#'
#' with(sleep,
#'      yuenTTest(extra[group == 1],
#'                extra[group == 2]))
#'
#' yuenTTest(extra ~ group, data = sleep)
#'
#' @rdname yuenTTest
#' @family test.location
#' @concept location-test
#' @concept robust-statistics
#'
#' @export
yuenTTest <- function(x, ...)
  UseMethod("yuenTTest")


# -------------------------------------------------------------------------
# Formula method
# -------------------------------------------------------------------------

#' @rdname yuenTTest
#' @export
yuenTTest.formula <- function(formula,
                              data,
                              subset,
                              na.action = na.pass,
                              ...) {
  
  if (missing(formula) || length(formula) != 3L)
    stop("'formula' missing or incorrect")
  
  args <- list(
    formula   = formula,
    na.action = na.action,
    allowed   = "two-sample-independent"
  )
  
  if (!missing(data))
    args$data <- data
  
  if (!missing(subset))
    args$subset <- substitute(subset)
  
  d <- do.call(
    resolveFormula,
    args
  )
  
  # resolveFormula() returns d$x as the FULL response (both groups,
  # length n) and d$group as the matching full-length factor - d$y is
  # only a convenience alias for group 2, never group 1. Split on
  # d$group explicitly instead of relying on d$x/d$y directly, or this
  # silently compares "all observations" against "group 2 only".
  groups <- split(d$x, d$group)
  
  res <- yuenTTest.default(
    x = groups[[1L]],
    y = groups[[2L]],
    ...
  )
  
  res$data.name <- d$data.name
  
  res
  
}


# -------------------------------------------------------------------------
# Default method
# -------------------------------------------------------------------------

#' @rdname yuenTTest
#' @export
yuenTTest.default <- function(
    x,
    y = NULL,
    alternative = c("two.sided", "less", "greater"),
    mu = 0,
    paired = FALSE,
    conf.level = 0.95,
    trim = 0.2,
    ...
) {
  
  alternative <- match.arg(alternative)
  
  if (!is.numeric(mu) ||
      length(mu) != 1L ||
      is.na(mu)) {
    stop("'mu' must be a single numeric value")
  }
  
  if (!is.numeric(conf.level) ||
      length(conf.level) != 1L ||
      !is.finite(conf.level) ||
      conf.level <= 0 ||
      conf.level >= 1) {
    stop("'conf.level' must be in (0,1)")
  }
  
  if (!is.numeric(trim) ||
      length(trim) != 1L ||
      is.na(trim) ||
      trim < 0 ||
      trim >= 0.5) {
    stop("'trim' must satisfy 0 <= trim < 0.5")
  }
  
  if (!is.null(y)) {
    
    dname <- paste(
      deparse1(substitute(x)),
      "and",
      deparse1(substitute(y))
    )
    
    if (paired) {
      
      ok <- complete.cases(x, y)
      
      x <- x[ok]
      y <- y[ok]
      
    } else {
      
      x <- x[is.finite(x)]
      y <- y[is.finite(y)]
    }
    
  } else {
    
    dname <- deparse1(substitute(x))
    
    if (paired)
      stop("'y' is missing for paired test")
    
    x <- x[is.finite(x)]
  }
  
  ## ---------------------------------------------------------------------
  ## One-sample / paired
  ## ---------------------------------------------------------------------
  
  if (is.null(y) || paired) {
    
    if (paired) {
      
      d <- x - y
      
      d <- d[is.finite(d)]
      
      n <- length(d)
      
      if (n < 2)
        stop("not enough paired observations")
      
      g <- floor(trim * n)
      
      df <- n - 2 * g - 1
      
      if (df <= 0)
        stop("trim level too large for sample size")
      
      md <- mean(d, trim = trim)
      
      se <- .trimmedSE(d, trim)
      
      if (se < 10 * .Machine$double.eps * abs(md))
        stop("data are essentially constant")
      
      tstat <- (md - mu) / se
      
      estimate <- c(
        "difference in trimmed means" = md
      )
      
      method <- "Yuen Paired-Sample Trimmed Mean t-test"
      
    } else {
      
      n <- length(x)
      
      if (n < 2)
        stop("not enough 'x' observations")
      
      g <- floor(trim * n)
      
      df <- n - 2 * g - 1
      
      if (df <= 0)
        stop("trim level too large for sample size")
      
      mx <- mean(x, trim = trim)
      
      se <- .trimmedSE(x, trim)
      
      if (se < 10 * .Machine$double.eps * abs(mx))
        stop("data are essentially constant")
      
      tstat <- (mx - mu) / se
      
      estimate <- c(
        "trimmed mean of x" = mx
      )
      
      method <- "Yuen One-Sample Trimmed Mean t-test"
    }
    
  } else {
    
    ## -------------------------------------------------------------------
    ## Two-sample
    ## -------------------------------------------------------------------
    
    nx <- length(x)
    ny <- length(y)
    
    if (nx < 2)
      stop("not enough 'x' observations")
    
    if (ny < 2)
      stop("not enough 'y' observations")
    
    gx <- floor(trim * nx)
    gy <- floor(trim * ny)
    
    dfx <- nx - 2 * gx - 1
    dfy <- ny - 2 * gy - 1
    
    if (dfx <= 0 || dfy <= 0)
      stop("trim level too large for sample size")
    
    mx <- mean(x, trim = trim)
    my <- mean(y, trim = trim)
    
    vx <- .winsorVar(x, trim)
    vy <- .winsorVar(y, trim)
    
    stderrx <- ((nx - 1) * vx) /
      ((dfx + 1) * dfx)
    
    stderry <- ((ny - 1) * vy) /
      ((dfy + 1) * dfy)
    
    se <- sqrt(stderrx + stderry)
    
    if (se < 10 * .Machine$double.eps *
        max(abs(mx), abs(my))) {
      stop("data are essentially constant")
    }
    
    df <- (stderrx + stderry)^2 /
      ((stderrx^2 / dfx) +
         (stderry^2 / dfy))
    
    tstat <- (mx - my - mu) / se
    
    estimate <- c(
      "trimmed mean of x" = mx,
      "trimmed mean of y" = my
    )
    
    method <- "Yuen Two-Sample Trimmed Mean t-test"
  }
  
  ## ---------------------------------------------------------------------
  ## Confidence interval
  ## ---------------------------------------------------------------------
  
  alpha <- 1 - conf.level
  
  if (alternative == "less") {
    
    pval <- stats::pt(tstat, df)
    
    cint <- c(
      -Inf,
      estimate[1] - mu +
        stats::qt(conf.level, df) * se
    )
    
  } else if (alternative == "greater") {
    
    pval <- stats::pt(
      tstat,
      df,
      lower.tail = FALSE
    )
    
    cint <- c(
      estimate[1] - mu -
        stats::qt(conf.level, df) * se,
      Inf
    )
    
  } else {
    
    pval <- 2 * stats::pt(
      -abs(tstat),
      df
    )
    
    crit <- stats::qt(
      1 - alpha / 2,
      df
    )
    
    cint <- estimate[1] - mu +
      c(-crit, crit) * se
  }
  
  names(cint) <- c("lower", "upper")
  
  attr(cint, "conf.level") <- conf.level
  
  ## ---------------------------------------------------------------------
  ## Return object
  ## ---------------------------------------------------------------------
  
  names(tstat) <- "t"
  
  rval <- list(
    
    statistic = tstat,
    
    parameter = c(
      df   = df,
      trim = trim
    ),
    
    p.value = pval,
    
    conf.int = cint,
    
    estimate = estimate,
    
    null.value = c(
      "trimmed mean difference" = mu
    ),
    
    alternative = alternative,
    
    method = method,
    
    data.name = dname
  )
  
  class(rval) <- "htest"
  
  rval
  
}


# == internal helper functions ===========================================


.trimmedSE <- function(z, trim) {
  
  q <- stats::quantile(
    z,
    probs = c(trim, 1 - trim),
    type = 7,
    na.rm = TRUE
  )
  
  wz <- winsorize(
    z,
    val = q
  )
  
  sqrt(stats::var(wz)) /
    ((1 - 2 * trim) * sqrt(length(z)))
}




.winsorVar <- function(z, trim) {
  
  q <- stats::quantile(
    z,
    probs = c(trim, 1 - trim),
    type = 7,
    na.rm = TRUE
  )
  
  wz <- winsorize(
    z,
    val = q
  )
  
  stats::var(wz)
}




