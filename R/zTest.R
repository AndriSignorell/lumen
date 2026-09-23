
#' Z-Test for Testing Means With Known Population Standard Deviations
#' 
#' A parametric test for the mean of a normal distribution when the population 
#' variance is known, or for comparing two means with known variances, based 
#' on the standard normal distribution.
#' 
#' Compute the test of hypothesis and compute confidence interval on the mean
#' of a population when the standard deviation of the population is known.
#' 
#' Most introductory statistical texts introduce inference by using the z-test
#' and z-based confidence intervals based on knowing the population standard
#' deviation. However statistical packages often do not include functions to do
#' z-tests since the t-test is usually more appropriate for real world
#' situations. This function is meant to be used during that short period of
#' learning when the student is learning about inference using z-procedures,
#' but has not learned the t-based procedures yet.  Once the student has
#' learned about the t-distribution the `t.test()` function should be used
#' instead of this one (but the syntax is very similar, so this function should
#' be an appropriate introductory step to learning `t.test()`).
#' 
#' The formula interface is only applicable for the 2-sample tests.
#' 
#' @name zTest
#' @aliases zTest zTest.default zTest.formula
#' @param x numeric vector of data values. Non-finite (e.g. infinite or
#' missing) values will be omitted.
#' @param y an optional numeric vector of data values: as with x non-finite
#' values will be omitted.
#' @param mu a number specifying the hypothesized mean of the population.
#' @param sd_pop a positive number specifying the known standard deviation
#' of the population. Required. For the two-sample test, this single value
#' is assumed to be the common known standard deviation of both
#' populations.
#' @param alternative a character string specifying the alternative hypothesis,
#' must be one of `"two.sided"` (default), `"greater"` or
#' `"less"`.  You can specify just the initial letter. \cr For one-sample
#' tests, `alternative` refers to the true mean of the parent population
#' in relation to the hypothesized value of the mean.
#' @param paired a logical indicating whether you want a paired z-test.
#' Only available in the default method: the formula interface describes
#' independent groups and does not identify pairs.
#' @param conf.level confidence level for the interval computation.
#' @param formula a formula of the form `lhs ~ rhs` where `lhs` gives
#' the data values and `rhs` a factor with two levels giving the
#' corresponding groups.
#' @param data an optional matrix or data frame (or similar: see
#' [model.frame()]) containing the variables in the formula
#' `formula`.  By default the variables are taken from
#' `environment(formula)`.
#' @param subset an optional vector specifying a subset of observations to be
#' used.
#' @param na.action a function which indicates what should happen when the data
#' contain `NA`s. Defaults to `getOption("na.action")`.
#' @param \dots further arguments to be passed to or from methods.
#' 
#' @return A list with class "`htest`" containing the following
#' components: \item{statistic}{ the value of the z-statistic.} \item{p.value}{
#' the p-value for the test} \item{conf.int}{a confidence interval for the mean
#' appropriate to the specified alternative hypothesis.} \item{estimate}{the
#' estimated mean or difference in means depending on whether it was a
#' one-sample test or a two-sample test.} \item{null.value}{the specified
#' hypothesized value of the mean or mean difference depending on whether it
#' was a one-sample test or a two-sample test.} \item{alternative}{a character
#' string describing the alternative hypothesis.} \item{method}{ a character
#' string indicating what type of test was performed.} \item{data.name}{a
#' character string giving the name(s) of the data.}
#' 
#' @seealso [t.test()], [print.htest()]
#' @references Stahel, W. (2002) *Statistische Datenanalyse, 4th ed*,
#' vieweg
#' 
#' @examples
#' 
#' x <- rnorm(25, 100, 5)
#' zTest(x, mu=99, sd_pop=5)
#' 
#' # the classic interface
#' with(sleep, zTest(extra[group==1], extra[group==2], sd_pop=2))
#' 
#' # the formula interface
#' zTest(extra ~ group, data=sleep, sd_pop=2)
#' 
#' 
#' # Stahel (2002), pp. 186, 196
#' 
#' Tyres <- data.frame(A=c(44.5,55,52.5,50.2,45.3,46.1,52.1,50.5,50.6,49.2),
#'                       B=c(44.9,54.8,55.6,55.2,55.6,47.7,53,49.1,52.3,50.7))
#' with(Tyres, zTest(A, B, sd_pop=3, paired=TRUE))
#' 
#' 
#' Oxen <- data.frame(ext=c(2.7,2.7,1.1,3.0,1.9,3.0,3.8,3.8,0.3,1.9,1.9),
#'                    int=c(6.5,5.4,8.1,3.5,0.5,3.8,6.8,4.9,9.5,6.2,4.1))
#' with(Oxen, zTest(int, ext, sd_pop=1.8, paired=FALSE))
#' 
#' @rdname zTest
#' @family test.location
#' @concept location-test
#' @concept parametric
#'
#' @export
zTest <- function (x, ...)
  UseMethod("zTest")



#' @rdname zTest
#' @export
zTest.formula <- function(formula,
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
  # caught here - passed on in '...', zTest.default() would match them
  # partially to its own 'paired' (a gap stats::t.test.formula() still has)
  if (!isFALSE(paired))
    stop("'paired' must be FALSE in the formula interface; ",
         "use zTest(x, y, paired = TRUE)")
  
  # direct call, never do.call(): do.call() evaluates the substituted
  # subset expression in this frame, where the data columns do not exist
  subset_expr <- if (!missing(subset)) substitute(subset) else NULL

  d <- resolveFormula(formula, data,
                      subset    = subset_expr,
                      na.action = na.action,
                      allowed   = "two-sample-independent")

  # d$x is the full response (both groups); d$y is only a convenience
  # alias for group 2. Split explicitly by d$group instead of relying
  # on d$x/d$y directly.
  groups <- split(d$x, d$group)

  res <- zTest.default(
    x = groups[[1L]],
    y = groups[[2L]],
    ...
  )
  
  res$data.name <- d$dataName
  res
}



#' @rdname zTest
#' @export
zTest.default <- function (x, y = NULL, alternative = c("two.sided", "less", "greater"),
                           paired = FALSE, mu = 0, sd_pop, conf.level = 0.95,  ...)  {
  
  alternative <- match.arg(alternative)

  if (!is.numeric(mu) || length(mu) != 1L || !is.finite(mu))
    stop("'mu' must be a single number (finite)")

  if (!is.logical(paired) || length(paired) != 1L || is.na(paired))
    stop("'paired' must be a single non-missing logical value")

  # the known standard deviation is the whole point of the z-test; checking
  # it here replaces the t.test-style "essentially constant" guard, which
  # makes no sense when the standard error does not depend on the data
  if (missing(sd_pop))
    stop("'sd_pop' (the known population standard deviation) is required")
  if (!is.numeric(sd_pop) || length(sd_pop) != 1L ||
      !is.finite(sd_pop) || sd_pop <= 0)
    stop("'sd_pop' must be a single positive number")

  # all-NA input (logical NA) is let through to the "not enough
  # observations" checks below
  .num <- function(v) is.numeric(v) || all(is.na(v))
  if (!.num(x) || (!is.null(y) && !.num(y)))
    stop("'x' and 'y' must be numeric")
  if (!missing(conf.level) && (length(conf.level) != 1 || !is.finite(conf.level) ||
                               conf.level < 0 || conf.level > 1))
    stop("'conf.level' must be a single number between 0 and 1")
  
  if (!is.null(y)) {
    dname <- paste(deparse1(substitute(x)), "and", deparse1(substitute(y)))
    
    if (paired) {
      if (length(x) != length(y))
        stop("'x' and 'y' must have the same length for a paired test")
      xok <- yok <- is.finite(x) & is.finite(y)
    } else {
      yok <- is.finite(y)
      xok <- is.finite(x)
    }
    
    y <- y[yok]
    
  } else {
    dname <- deparse1(substitute(x))
    if (paired)
      stop("'y' is missing for paired test")
    xok <- is.finite(x)
    yok <- NULL
  }
  x <- x[xok]
  
  if (paired) {
    x <- x - y
    y <- NULL
  }
  
  nx <- length(x)
  mx <- mean(x)
  
  if (is.null(y)) {
    # with a known sd_pop a single observation suffices; n >= 2 was
    # inherited from t.test(), which has to estimate the standard deviation
    if (nx < 1)
      stop("not enough 'x' observations")
    stderr <- sd_pop / sqrt(nx)
    zstat <- (mx - mu)/stderr
    
    method <- if (paired)
      "Paired z-test" else "One Sample z-test"
    estimate <- setNamesX(mx, if (paired)
      "mean of the differences"
      else "mean of x")
  }
  else {
    ny <- length(y)
    if (nx < 1)
      stop("not enough 'x' observations")
    if (ny < 1)
      stop("not enough 'y' observations")
    my <- mean(y)
    
    method <- paste("Two Sample z-test")
    estimate <- c(mx, my)
    names(estimate) <- c("mean of x", "mean of y")
    
    stderr <- sd_pop * sqrt(1/nx + 1/ny)
    zstat <- (mx - my - mu)/stderr
  }
  if (alternative == "less") {
    pval <- pnorm(zstat)
    cint <- c(-Inf, zstat + qnorm(conf.level))
  }
  else if (alternative == "greater") {
    pval <- pnorm(zstat, lower.tail = FALSE)
    cint <- c(zstat - qnorm(conf.level), Inf)
  }
  else {
    pval <- 2 * pnorm(-abs(zstat))
    alpha <- 1 - conf.level
    cint <- qnorm(1 - alpha/2)
    cint <- zstat + c(-cint, cint)
  }
  cint <- mu + cint * stderr
  names(zstat) <- "z"
  names(mu) <- if (paired || !is.null(y))
    "difference in means"
  else "mean"
  names(sd_pop) <- "Std. Dev. Population"
  attr(cint, "conf.level") <- conf.level
  rval <- list(
    statistic = zstat, parameter = sd_pop, p.value = pval,
    conf.int = cint, estimate = estimate, null.value = mu, stderr = stderr,
    alternative = alternative, method = method, data.name = dname )
  class(rval) <- "htest"
  return(rval)
}


