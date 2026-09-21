
#' Simple Bootstrap Confidence Intervals 
#' 
#' Convenience wrapper for calculating bootstrap confidence intervals for
#' univariate and bivariate statistics. 
#' 
#' 
#' @param x a (non-empty) numeric vector of data values.
#' @param y NULL (default) or a vector with compatible dimensions to `x`,
#' when a bivariate statistic is used.
#' @param FUN the function to be used.
#' @param bci.method a vector of character strings representing the type of
#' intervals required. The value should be any subset of the values
#' `"norm"`, `"basic"`, `"stud"`, `"perc"`, `"bca"`,
#' as it is passed on as `method` to [boot::boot.ci()].
#' @param conf.level confidence level of the interval.
#' @param sides a character string specifying the side of the confidence
#' interval, must be one of `"two.sided"` (default), `"left"` or
#' `"right"`. You can specify just the initial letter. `"left"` would
#' be analogue to a hypothesis of `"greater"` in a `t.test`.
#' @param ... further arguments are passed to the function `FUN`.
#' @param R number of bootstrap replicates. Usually this will be a single
#' positive integer. For importance resampling, some resamples may use one set
#' of weights and others use a different set of weights. In this case `R`
#' would be a vector of integers where each component gives the number of
#' resamples from each of the rows of weights.
#' 
#' @return A named numeric vector with three elements:
#' \describe{
#'   \item{`est`}{the estimate calculated by `FUN`.}
#'   \item{`lci`}{lower confidence interval bound.}
#'   \item{`uci`}{upper confidence interval bound.}
#' }
#' 
#' @examples
#' 
#' set.seed(1984)
#' bootCI(mtcars$mpg, FUN=mean, na.rm=TRUE, bci.method="basic")
#' bootCI(mtcars$mpg, FUN=mean, trim=0.1, na.rm=TRUE, bci.method="basic")
#' 
#' # bootCI(mtcars$mpg, FUN=DescToolsX::skewX, na.rm=TRUE, bci.method="basic")
#' 
#' # bootCI(Pizza$operator, Pizza$area, FUN=cramerV)
#' 
#' spearman <- function(x,y) cor(x, y, method="spearman", use="p")
#' bootCI(mtcars$mpg, mtcars$hp, FUN=spearman)
#' 
#' 
#' 
#' @family ci.general  
#' @concept confidence-interval  
#' @concept bootstrap
#'
#'
#' @export
bootCI <- function(x, y=NULL, FUN, ..., bci.method = c("norm", "basic", "stud", "perc", "bca"),
                   conf.level = 0.95, sides = c("two.sided", "left", "right"), R = 999) {

  # evaluated here, in the caller's frame: substitute() handed unevaluated
  # expressions to do.call(), which then resolved them inside the boot()
  # statistic, where no caller variable is visible
  dots <- list(...)
  bci.method <- match.arg(bci.method)
  sides <- match.arg(sides)

  if (sides != "two.sided") {
    if (conf.level <= 0.5)
      stop(gettextf("a one-sided interval needs 'conf.level' above 0.5, not %g",
                    conf.level), domain = NA)
    conf.level <- 1 - 2 * (1 - conf.level)
  }

  stat <- if (!is.null(y)) {
    function(x, d) do.call(FUN, c(list(x[d], y[d]), dots))
  } else if (is.matrix(x) || is.data.frame(x)) {
    function(x, d) do.call(FUN, c(list(x[d, , drop = FALSE]), dots))
  } else {
    function(x, d) do.call(FUN, c(list(x[d]), dots))
  }

  boot.fun <- boot::boot(x, stat, R = R)

  ci <- boot::boot.ci(boot.fun, conf = conf.level, type = bci.method)

  # by name, not ci[[4]]: a dropped component ('stud' without variances)
  # made the positional access fail with "subscript out of bounds"
  bnd <- .bootCIBounds(ci, bci.method)

  res <- c(est = unname(boot.fun$t0[1L]), lci = bnd[1L], uci = bnd[2L])

  if (sides == "left")
    res[["uci"]] <- Inf
  else if (sides == "right")
    res[["lci"]] <- -Inf

  res
}
