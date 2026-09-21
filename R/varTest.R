
#' Variance Test for Testing One Variance or Comparing Two Variances
#'
#' Performs a one-sample or two-sample test for variance, analogous to
#' [t.test()], with support for classical and likelihood-based
#' lowest-density (LD) two-sided p-values.
#'
#' @param x a numeric vector of data values, or a formula.
#' @param y an optional second numeric vector. If provided, a two-sample variance
#' test is performed.
#' @param sigma2_0 a numeric value specifying the null hypothesis variance for
#' the one-sample test. Required if `y` is `NULL`.
#' @param alternative character string specifying the alternative hypothesis.
#' Must be one of `"two.sided"`, `"less"`, or `"greater"`.
#' @param type character string specifying the test type:
#' \itemize{
#'   \item `"classic"`: uses the conventional two-sided p-value
#'   \eqn{2 \cdot \min(P(T \le t), P(T \ge t))}.
#'   \item `"ld"`: uses a likelihood-based lowest-density (LD) definition,
#'   i.e. the probability of observing values with density less than or equal
#'   to the observed density under the null distribution.
#' }
#' @param ... further arguments passed to methods.
#'
#' @name varTest
#' @details
#' The null hypothesis is that the ratio of the variances of the populations
#' from which `x` and `y` were drawn, or in the data to which the
#' linear models `x` and `y` were fitted, is equal to `ratio`.
#' 
#' For the one-sample test, the test statistic follows a chi-squared distribution:
#' \deqn{X^2 = (n - 1) S^2 / \sigma_0^2}
#'
#' For the two-sample test, the test statistic follows an F distribution:
#' \deqn{F = S_x^2 / S_y^2}
#'
#' The LD method corresponds to a likelihood-ratio interpretation of extremeness,
#' which is particularly appropriate for asymmetric null distributions such as
#' chi-squared and F distributions.
#'
#' The formula interface is only applicable for the 2-sample tests.
#' 
#' @inheritParams Formulas
#' 
#' @return An object of class `"htest"` with components:
#' \item{statistic}{the test statistic.}
#' \item{parameter}{degrees of freedom.}
#' \item{p.value}{the p-value of the test.}
#' \item{estimate}{sample variance(s).}
#' \item{null.value}{the null hypothesis value (one-sample only).}
#' \item{alternative}{the alternative hypothesis.}
#' \item{method}{a character string indicating the test performed.}
#' \item{data.name}{description of the data.}
#'
#' @seealso [var.test()], [bartlett.test()] for testing
#' homogeneity of variances in more than two samples from normal distributions;
#' [ansari.test()] and [mood.test()] for two rank based
#' (nonparametric) two-sample tests for difference in scale.
#' 
#' @examples
#' set.seed(1)
#' x <- rnorm(20, sd = 3)
#' y <- rnorm(25, sd = 2)
#'
#' # One-sample test
#' varTest(x, sigma2_0 = 9, type = "classic")
#' varTest(x, sigma2_0 = 9, type = "ld")
#'
#' # Two-sample test
#' varTest(x, y, type = "classic")
#' varTest(x, y, type = "ld")
#'
#' # Formula interface
#' df <- data.frame(
#'   value = c(x, y),
#'   group = rep(c("A", "B"), c(length(x), length(y)))
#' )
#' varTest(value ~ group, data = df)
#'
#' @rdname varTest
#' @family test.variance
#' @concept variance-test
#'
#' @export
varTest <- function(x, ...) UseMethod("varTest")


#' @rdname varTest
#' @export
varTest.default <- function(x, y = NULL, sigma2_0 = NULL,
                            alternative = c("two.sided", "less", "greater"),
                            type = c("classic", "ld"), ...) {
  
  
  # https://stats.stackexchange.com/questions/140107/p-value-in-a-two-tail-test-with-asymmetric-null-distribution
  
  # https://stats.stackexchange.com/questions/195469/calculating-p-values-for-two-tail-test-for-population-variance
  
  # What you are dealing with in this question is a two-sided 
  # variance test, which is a specific case of a two-sided test 
  # with an asymmetric null distribution. The p-value is the total 
  # area under the null density for all values in the lower and 
  # upper tails of that density that are at least as "extreme" 
  # (i.e., at least as conducive to the alternative hypothesis) 
  # as the observed test statistic. Because this test has an 
  # asymmetric null distribution, we need to specify exactly 
  # what we mean by "extreme".
  # 
  # Lowest-density p-value calculation: The most sensible thing 
  # method of two-sided hypothesis testing is to interpret 
  # "more extreme" as meaning a lower value of the null density. 
  # This is the interpretation used in a standard likelihood-ratio 
  # (LR) test. Under this method , the p-value is the probability of 
  # falling in the "lowest density region", where the density 
  # cut-off is the density at the observed test statistic. With 
  # an asymmetric null distribution, this leads you to a p-value 
  # calculated with unequal tails.
  
  
  alternative <- match.arg(alternative)
  type        <- match.arg(type)

  if (!is.numeric(x) || length(x) < 2L)
    stop("'x' must be a numeric vector with at least two observations")

  # =============================
  # One-sample test
  # =============================
  if (is.null(y)) {

    if (is.null(sigma2_0))
      stop("sigma2_0 must be provided for one-sample test.")

    if (!is.numeric(sigma2_0) || length(sigma2_0) != 1L ||
        !is.finite(sigma2_0) || sigma2_0 <= 0)
      stop("'sigma2_0' must be a single positive finite number")

    nu   <- length(x) - 1
    s2   <- var(x)
    stat <- c("X-squared" = nu * s2 / sigma2_0)

    pdist <- function(q, lower.tail = TRUE) pchisq(q, nu, lower.tail = lower.tail)
    ddist <- function(q) dchisq(q, nu)
    mode  <- max(nu - 2, 0)

    parameter  <- c(df = nu)
    estimate   <- c(variance = s2)
    null.value <- c(variance = sigma2_0)
    method     <- paste0("One-sample variance test (", type, ")")
    data.name  <- deparse1(substitute(x))

    # =============================
    # Two-sample test
    # =============================
  } else {

    if (!is.numeric(y) || length(y) < 2L)
      stop("'y' must be a numeric vector with at least two observations")

    nu1  <- length(x) - 1
    nu2  <- length(y) - 1
    s2x  <- var(x)
    s2y  <- var(y)
    stat <- c(F = s2x / s2y)

    pdist <- function(q, lower.tail = TRUE) pf(q, nu1, nu2, lower.tail = lower.tail)
    ddist <- function(q) df(q, nu1, nu2)
    mode  <- if (nu1 > 2) (nu1 - 2) / nu1 * nu2 / (nu2 + 2) else 0

    parameter  <- c(df1 = nu1, df2 = nu2)
    estimate   <- c("var(x)" = s2x, "var(y)" = s2y)
    null.value <- NULL
    method     <- paste0("Two-sample variance test (", type, ")")
    data.name  <- paste(deparse1(substitute(x)), "and", deparse1(substitute(y)))
  }

  # Filter() drops null.value for the two-sample test
  structure(
    Filter(Negate(is.null), list(
         statistic   = stat,
         parameter   = parameter,
         p.value     = .varTestPValue(unname(stat), pdist, ddist, mode,
                                      alternative, type),
         estimate    = estimate,
         null.value  = null.value,
         alternative = alternative,
         method      = method,
         data.name   = data.name)),
    class = "htest")
}



#' @rdname varTest
#' @export
varTest.formula <- function(formula,
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

  # quote = TRUE: without it do.call() evaluates the subset expression in
  # this frame instead of passing it on unevaluated, so a subset naming a
  # column of 'data' fails with "object not found"
  d <- do.call(resolveFormula, args, quote = TRUE)

  if (nlevels(d$group) != 2L)
    stop("grouping factor must have exactly 2 levels")

  groups <- split(d$x, d$group)

  res <- varTest.default(
    x = groups[[1L]],
    y = groups[[2L]],
    ...
  )

  res$data.name <- d$dataName      # htest uses data.name, resolveFormula dataName
  res
}



# == internal helper functions ================================================


.varTestPValue <- function(q, pdist, ddist, mode, alternative, type) {

  if (is.na(q))
    return(NA_real_)

  switch(alternative,
         less      = pdist(q),
         greater   = pdist(q, lower.tail = FALSE),
         two.sided = if (type == "classic")
           2 * min(pdist(q), pdist(q, lower.tail = FALSE))
         else
           .ldPValue(q, pdist, ddist, mode))
}


.ldPValue <- function(q, pdist, ddist, mode, tol = 1e-12) {

  # Lowest-density two-sided p-value P(d(T) <= d(q)) for a unimodal null
  # density d with interior mode 'mode'. mode = 0 means d is decreasing on
  # (0, Inf) (chi-squared df <= 2, F df1 <= 2): no second tail, the LD
  # region is the plain upper tail.
  if (mode <= 0)
    return(pdist(q, lower.tail = FALSE))

  if (q == mode)
    return(1)

  dObs <- ddist(q)
  g    <- function(t) ddist(t) - dObs

  # d(0) = 0 whenever mode > 0, so [0, mode] always brackets the left root.
  # Right of the mode d is decreasing, so the upper end is extended until
  # the sign changes. (The former fixed brackets [1e-10, mode] and
  # [mode, q(0.999999)] failed for extreme statistics, e.g. chi-squared
  # df = 3, q = 30.) A tight tol matters: the default ~1.2e-4 cost up to
  # 4e-6 in the p-value.
  if (q > mode) {
    root <- uniroot(g, c(0, mode), tol = tol)$root
    pdist(root) + pdist(q, lower.tail = FALSE)
  } else {
    root <- uniroot(g, c(mode, 2 * mode + 1), extendInt = "downX",
                    tol = tol)$root
    pdist(q) + pdist(root, lower.tail = FALSE)
  }
}
