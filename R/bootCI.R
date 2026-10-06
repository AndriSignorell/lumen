
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
#' @param conf.level confidence level of the interval.
#' @param sides a character string specifying the side of the confidence
#' interval, must be one of `"two.sided"` (default), `"left"` or
#' `"right"`. You can specify just the initial letter. `"left"` would
#' be analogue to a hypothesis of `"greater"` in a `t.test`.
#' @param R number of bootstrap replicates, a single positive whole number.
#' @param ... further arguments. The bootstrap options are taken out first,
#' as in the other interval functions of the package: `type`, the interval
#' type passed to [boot::boot.ci()], one of `"bca"` (default), `"perc"`,
#' `"basic"`, `"norm"` or `"stud"`, and `parallel` and `ncpus`, passed to
#' [boot::boot()]. Everything else is passed to `FUN`.
#'
#' @details
#' `type`, `parallel` and `ncpus` therefore cannot reach `FUN` through the
#' dots. A statistic that has an argument of one of these names - the
#' `type` of [quantile()], say - is wrapped:
#' `FUN = function(z) quantile(z, 0.9, type = 6)`.
#'
#' `"stud"` needs a variance estimate for every replicate, which a general
#' `FUN` does not deliver; [boot::boot.ci()] then returns no such interval
#' and `bootCI()` stops with a message saying so.
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
#' bootCI(mtcars$mpg, FUN=mean, na.rm=TRUE)
#' bootCI(mtcars$mpg, FUN=mean, trim=0.1, na.rm=TRUE, type="basic")
#' 
#' # bootCI(mtcars$mpg, FUN=DescToolsX::skewX, na.rm=TRUE, type="basic")
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
bootCI <- function(x, y=NULL, FUN, conf.level = 0.95,
                   sides = c("two.sided", "left", "right"), R = 999, ...) {

  # evaluated here, in the caller's frame: substitute() handed unevaluated
  # expressions to do.call(), which then resolved them inside the boot()
  # statistic, where no caller variable is visible
  dots <- list(...)
  sides <- match.arg(sides)
  checkConfLevel(conf.level, allowNA = FALSE)

  # the bootstrap options leave the dots, the rest belongs to FUN. Split by
  # a logical index: setdiff() on the names would drop unnamed elements.
  nms <- names(dots)
  if (is.null(nms)) nms <- rep("", length(dots))
  isBoot   <- nms %in% .bootArgNames
  bootArgs <- .extractBootArgs(c(list(R = R), dots[isBoot]))
  dots     <- dots[!isBoot]

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

  boot.fun <- boot::boot(x, stat, R = bootArgs$R,
                         parallel = bootArgs$parallel, ncpus = bootArgs$ncpus)

  ci <- boot::boot.ci(boot.fun, conf = conf.level, type = bootArgs$type)

  # by name, not ci[[4]]: a dropped component ('stud' without variances)
  # made the positional access fail with "subscript out of bounds"
  bnd <- .bootCIBounds(ci, bootArgs$type)

  res <- c(est = unname(boot.fun$t0[1L]), lci = bnd[1L], uci = bnd[2L])

  if (sides == "left")
    res[["uci"]] <- Inf
  else if (sides == "right")
    res[["lci"]] <- -Inf

  res
}
