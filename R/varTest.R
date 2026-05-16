
#' Variance Test (One- and Two-Sample) with Classic and LD Methods
#'
#' Performs a one-sample or two-sample test for variance, analogous to
#' \code{\link[stats]{t.test}}, with support for classical and likelihood-based
#' lowest-density (LD) two-sided p-values.
#'
#' @param x A numeric vector of data values, or a formula.
#' @param y An optional second numeric vector. If provided, a two-sample variance
#' test is performed.
#' @param sigma2_0 A numeric value specifying the null hypothesis variance for
#' the one-sample test. Required if \code{y} is \code{NULL}.
#' @param alternative Character string specifying the alternative hypothesis.
#' Must be one of \code{"two.sided"}, \code{"less"}, or \code{"greater"}.
#' @param type Character string specifying the test type:
#' \itemize{
#'   \item \code{"classic"}: uses the conventional two-sided p-value
#'   \eqn{2 \cdot \min(P(T \le t), P(T \ge t))}.
#'   \item \code{"ld"}: uses a likelihood-based lowest-density (LD) definition,
#'   i.e. the probability of observing values with density less than or equal
#'   to the observed density under the null distribution.
#' }
#' @param ... Further arguments passed to methods.
#'
#' @name varTest
#' @details
#' The null hypothesis is that the ratio of the variances of the populations
#' from which \code{x} and \code{y} were drawn, or in the data to which the
#' linear models \code{x} and \code{y} were fitted, is equal to \code{ratio}.
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
#' @return An object of class \code{"htest"} with components:
#' \item{statistic}{The test statistic.}
#' \item{parameter}{Degrees of freedom.}
#' \item{p.value}{The p-value of the test.}
#' \item{estimate}{Sample variance(s).}
#' \item{null.value}{The null hypothesis value (one-sample only).}
#' \item{alternative}{The alternative hypothesis.}
#' \item{method}{A character string indicating the test performed.}
#' \item{data.name}{Description of the data.}
#'
#' @seealso \code{\link{var.test}}, \code{\link{bartlett.test}} for testing
#' homogeneity of variances in more than two samples from normal distributions;
#' \code{\link{ansari.test}} and \code{\link{mood.test}} for two rank based
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
#' @concept hypothesis-testing
#' @concept descriptive-statistics
#' @concept confidence-intervals
#'
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
  
  # example:
  #   TwoSided_pval(15.35667, DF=17)
  #   TwoSided_pval(14.6489, DF=17)
  
  
  alternative <- match.arg(alternative)
  type <- match.arg(type)
  

  # =============================
  # One-sample test
  # =============================
  if (is.null(y)) {
    
    if (is.null(sigma2_0)) {
      stop("sigma2_0 must be provided for one-sample test.")
    }
    
    n <- length(x)
    df <- n - 1
    s2 <- var(x)
    x_obs <- df * s2 / sigma2_0
    
    if (type == "classic") {
      if (alternative == "two.sided") {
        p <- 2 * min(pchisq(x_obs, df), 1 - pchisq(x_obs, df))
      } else if (alternative == "less") {
        p <- pchisq(x_obs, df)
      } else {
        p <- 1 - pchisq(x_obs, df)
      }
    } else {
      if (alternative == "two.sided") {
        p <- .pchisqLD(x_obs, df)
      } else if (alternative == "less") {
        p <- pchisq(x_obs, df)
      } else {
        p <- 1 - pchisq(x_obs, df)
      }
    }
    
    result <- list(
      statistic = c("X-squared" = x_obs),
      parameter = c("df" = df),
      p.value = p,
      estimate = c("variance" = s2),
      null.value = c("variance" = sigma2_0),
      alternative = alternative,
      method = paste("One-sample variance test (", type, ")", sep = ""),
      data.name = deparse(substitute(x))
    )
    
    # =============================
    # Two-sample test
    # =============================
  } else {
    
    nx <- length(x)
    ny <- length(y)
    
    df1 <- nx - 1
    df2 <- ny - 1
    
    s2x <- var(x)
    s2y <- var(y)
    
    f_obs <- s2x / s2y
    
    if (type == "classic") {
      if (alternative == "two.sided") {
        p <- 2 * min(pf(f_obs, df1, df2), 1 - pf(f_obs, df1, df2))
      } else if (alternative == "less") {
        p <- pf(f_obs, df1, df2)
      } else {
        p <- 1 - pf(f_obs, df1, df2)
      }
    } else {
      if (alternative == "two.sided") {
        p <- .pfLD(f_obs, df1, df2)
      } else if (alternative == "less") {
        p <- pf(f_obs, df1, df2)
      } else {
        p <- 1 - pf(f_obs, df1, df2)
      }
    }
    
    result <- list(
      statistic = c("F" = f_obs),
      parameter = c("df1" = df1, "df2" = df2),
      p.value = p,
      estimate = c("var(x)" = s2x, "var(y)" = s2y),
      alternative = alternative,
      method = paste("Two-sample variance test (", type, ")", sep = ""),
      data.name = paste(deparse(substitute(x)), "and", deparse(substitute(y)))
    )
  }
  
  class(result) <- "htest"
  return(result)
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
    allowed   = "two.sample.independent"
  )
  
  if (!missing(data))
    args$data <- data
  
  if (!missing(subset))
    args$subset <- substitute(subset)
  
  d <- do.call(bedrock::resolveFormula, args)
  
  if (nlevels(d$group) != 2L)
    stop("grouping factor must have exactly 2 levels")
  
  groups <- split(d$x, d$group)
  
  res <- varTest.default(
    x = groups[[1L]],
    y = groups[[2L]],
    ...
  )
  
  res$data.name <- d$data.name
  res
}



# == internal helper functions ================================================


.pchisqLD <- function(x_obs, df) {
  
  # LD p-value chisq
  
  d_obs <- dchisq(x_obs, df)
  mode <- if (df >= 2) df - 2 else 0
  f <- function(x) dchisq(x, df) - d_obs
  
  if (x_obs >= mode) {
    root <- uniroot(f, c(1e-10, x_obs))$root
    pchisq(root, df) + (1 - pchisq(x_obs, df))
  } else {
    ub <- qchisq(0.999999, df)
    root <- uniroot(f, c(x_obs, ub))$root
    pchisq(x_obs, df) + (1 - pchisq(root, df))
  }
}


.pfLD <- function(f_obs, df1, df2) {
  
  # LD p-value F
  
  d_obs <- df(f_obs, df1, df2)
  
  mode <- if (df1 > 2) {
    (df1 - 2)/df1 * (df2/(df2 + 2))
  } else 0
  
  f_fun <- function(x) df(x, df1, df2) - d_obs
  
  if (f_obs >= mode) {
    root <- uniroot(f_fun, c(1e-10, f_obs))$root
    pf(root, df1, df2) + (1 - pf(f_obs, df1, df2))
  } else {
    ub <- qf(0.999999, df1, df2)
    root <- uniroot(f_fun, c(f_obs, ub))$root
    pf(f_obs, df1, df2) + (1 - pf(root, df1, df2))
  }
}




