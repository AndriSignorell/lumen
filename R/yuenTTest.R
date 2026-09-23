
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
#' For paired tests, trimming is performed on the paired differences, i.e.
#' the one-sample trimmed t-test (Tukey & McLaughlin, 1963) is applied to
#' \eqn{x - y}. This tests the trimmed mean of the differences, which in
#' general is not the difference of the trimmed means; the latter is what
#' e.g. `WRS2::yuend()` compares.
#'
#' @details
#' **Winsorizing.** With \eqn{g = \lfloor \mathrm{trim} \cdot n \rfloor},
#' the \eqn{g} smallest observations are set to the \eqn{(g+1)}-th order
#' statistic and the \eqn{g} largest to the \eqn{(n-g)}-th, as in Yuen
#' (1974) and Wilcox (2005). The winsorized variance must be taken this way,
#' at order statistics rather than at interpolated quantiles: the standard
#' error \eqn{\sqrt{(n-1) s_w^2 / (h(h-1))}} with \eqn{h = n - 2g} is
#' derived for exactly \eqn{g} replaced values in each tail, the same
#' \eqn{g} observations that [mean()] with `trim` removes. The results agree
#' with `WRS2::yuen()` and `PairedData::yuen.t.test()`.
#'
#' **Standard error.** In all three designs the squared standard error of a
#' trimmed mean is \eqn{(n-1) s_w^2 / (h(h-1))}, with \eqn{s_w^2} the
#' winsorized variance and \eqn{h} the number of observations left after
#' trimming; the degrees of freedom are \eqn{h - 1} (combined by Welch's
#' formula in the two-sample case). With `trim = 0`, or whenever
#' \eqn{g = 0}, the three tests reduce exactly to the corresponding
#' [t.test()]: one-sample, paired, and Welch. The one-sample version in
#' `WRS2::trimse()` uses the asymptotically equivalent
#' \eqn{s_w / ((1 - 2\,\mathrm{trim}) \sqrt{n})}; it differs in small
#' samples, where it inflates the standard error even if no observation is
#' trimmed (e.g. \eqn{n = 4}, `trim = 0.2`).
#'
#' The confidence interval is for the estimated parameter itself (the
#' trimmed mean, or the difference of trimmed means), independent of `mu`.
#'
#' @name yuenTTest
#' @aliases yuenTTest yuenTTest.default yuenTTest.formula
#'
#' @param x numeric vector of observations. Non-finite values (`NA`,
#'   `NaN`, `Inf`, `-Inf`) are removed; in the paired case the pair is
#'   removed.
#' @param y optional second numeric vector.
#' @param alternative character string specifying the alternative
#'   hypothesis. One of `"two.sided"`, `"less"`,
#'   or `"greater"`.
#' @param paired logical indicating whether a paired test is performed.
#'   Only available in the default method: the formula interface describes
#'   independent groups and does not identify pairs.
#' @param mu hypothesized trimmed mean (or trimmed mean difference).
#' @param conf.level confidence level for the confidence interval.
#' @param trim fraction of observations trimmed from each tail.
#'   Must satisfy `0 <= trim < 0.5`.
#' @param formula a formula of the form `lhs ~ rhs`.
#' @param data optional data frame for the formula interface.
#' @param subset optional subset expression.
#' @param na.action NA handling function.
#' @param \dots further arguments passed to methods.
#'
#' @return
#' An object of class `"htest"`.
#'
#' @seealso [t.test()]
#'
#' @references
#' Wilcox, R. R. (2005).
#' *Introduction to Robust Estimation and Hypothesis Testing*.
#' Academic Press.
#'
#' Tukey, J. W., & McLaughlin, D. H. (1963).
#' Less vulnerable confidence and significance procedures for location based
#' on a single sample: trimming/winsorization 1.
#' *Sankhya A*, 25, 331--352.
#'
#' Yuen, K. K. (1974).
#' The two-sample trimmed t for unequal population variances.
#' *Biometrika*, 61, 165--170.
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
                              paired = FALSE,
                              ...) {
  
  if (missing(formula) || length(formula) != 3L)
    stop("'formula' missing or incorrect")

  # the groups are split from independent rows; pairing them by position
  # would make the result depend on the row order within each group.
  # 'paired' is a formal argument so that abbreviations (pair = TRUE) are
  # caught here - passed on in '...', yuenTTest.default() would match them
  # partially to its own 'paired' (a gap stats::t.test.formula() still has)
  if (!isFALSE(paired))
    stop("'paired' must be FALSE in the formula interface; ",
         "use yuenTTest(x, y, paired = TRUE)")
  
  args <- list(
    formula   = formula,
    na.action = na.action,
    allowed   = "two-sample-independent"
  )
  
  if (!missing(data))
    args$data <- data
  
  if (!missing(subset))
    args$subset <- substitute(subset)
  
  d <- do.call(resolveFormula, args, quote = TRUE)
  
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
  
  res$data.name <- d$dataName
  
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
      !is.finite(mu)) {
    stop("'mu' must be a single finite numeric value")
  }

  if (!is.logical(paired) || length(paired) != 1L || is.na(paired))
    stop("'paired' must be a single non-missing logical value")
  
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
      
      if (se <= 10 * .Machine$double.eps * abs(md))
        stop("data are essentially constant")
      
      tstat <- (md - mu) / se
      
      est <- md

      estimate <- c(
        "trimmed mean of the differences" = md
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
      
      if (se <= 10 * .Machine$double.eps * abs(mx))
        stop("data are essentially constant")
      
      tstat <- (mx - mu) / se
      
      est <- mx

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
    
    # squared standard errors of the two trimmed means
    stderrx <- .trimmedSE(x, trim)^2
    stderry <- .trimmedSE(y, trim)^2
    
    se <- sqrt(stderrx + stderry)
    
    # <= rather than <: with se = 0 and trimmed means of 0 the bound is 0
    if (se <= 10 * .Machine$double.eps *
        max(abs(mx), abs(my))) {
      stop("data are essentially constant")
    }
    
    df <- (stderrx + stderry)^2 /
      ((stderrx^2 / dfx) +
         (stderry^2 / dfy))
    
    est <- mx - my

    tstat <- (est - mu) / se
    
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
      est + stats::qt(conf.level, df) * se
    )
    
  } else if (alternative == "greater") {
    
    pval <- stats::pt(
      tstat,
      df,
      lower.tail = FALSE
    )
    
    cint <- c(
      est - stats::qt(conf.level, df) * se,
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
    
    cint <- est + c(-crit, crit) * se
  }
  
  # the interval is for est (trimmed mean, trimmed mean of differences or
  # difference of trimmed means) - not for estimate[1], which in the
  # two-sample case is the trimmed mean of x alone, and not shifted by mu
  cint <- unname(cint)
  names(cint) <- c("lower", "upper")
  
  attr(cint, "conf.level") <- conf.level
  
  ## ---------------------------------------------------------------------
  ## Return object
  ## ---------------------------------------------------------------------
  
  names(tstat) <- "t"

  nullValue <- mu
  names(nullValue) <- if (is.null(y)) "trimmed mean" else "trimmed mean difference"
  
  rval <- list(
    
    statistic = tstat,
    
    parameter = c(
      df   = df,
      trim = trim
    ),
    
    p.value = pval,
    
    conf.int = cint,
    
    estimate = estimate,
    
    null.value = nullValue,
    
    alternative = alternative,
    
    method = method,
    
    data.name = dname
  )
  
  class(rval) <- "htest"
  
  rval
  
}


# == internal helper functions ===========================================


# winsorized variance as in Yuen (1974) / Wilcox (2005): the g = floor(trim*n)
# smallest values are set to the (g+1)-th order statistic, the g largest to
# the (n-g)-th. Interpolated quantiles (quantile(type = 7)) would winsorize
# at points between order statistics and no longer match h = n - 2g in the
# standard error, nor the g observations that mean(x, trim) removes.
.winsorVar <- function(z, trim) {

  z <- sort(z[!is.na(z)])
  n <- length(z)
  g <- floor(trim * n)

  if (g > 0L) {
    z[seq_len(g)]         <- z[g + 1L]
    z[(n - g + 1L):n]     <- z[n - g]
  }

  stats::var(z)
}


# standard error of the trimmed mean, sqrt((n-1) s_w^2 / (h (h-1))) with
# h = n - 2g the number of observations kept - the same formula in all three
# designs, and exactly the t-test standard error when g = 0. WRS2::trimse()
# uses s_w / ((1 - 2 trim) sqrt(n)) instead, which inflates the SE by
# 1 / (1 - 2 trim) even when nothing is trimmed (n = 4, trim = 0.2).
# Callers ensure h > 1.
.trimmedSE <- function(z, trim) {
  n <- sum(!is.na(z))
  h <- n - 2 * floor(trim * n)
  sqrt((n - 1) * .winsorVar(z, trim) / (h * (h - 1)))
}
