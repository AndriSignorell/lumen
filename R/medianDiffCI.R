
#' Confidence Interval for the Difference of Two Medians
#'
#' Calculates the confidence interval for the difference of two medians,
#' `median(x) - median(y)`, either with the distribution-free method of
#' Bonett and Price or with the bootstrap.
#'
#' `method = "bonett-price"` (default) estimates the standard error of
#' each sample median from a pair of order statistics,
#' \deqn{\widehat{se}_j = \frac{Y_{(n_j - c_j + 1)} - Y_{(c_j)}}{2 z_j},
#'   \quad c_j = \mathrm{round}\left(\frac{n_j}{2} - \sqrt{n_j}\right),}
#' where \eqn{z_j} is the normal quantile matching the exact binomial
#' coverage of \eqn{[Y_{(c_j)}, Y_{(n_j - c_j + 1)}]}, and combines both
#' into a normal interval. No distributional assumption is made, and the
#' variances of the two groups may differ. It is available for
#' independent samples only.
#'
#' `method = "boot"` resamples within each group (stratified), or whole
#' pairs if `paired = TRUE`. The bootstrap distribution of a median is
#' discrete and lumpy in small samples; `type = "perc"` is the robust
#' choice there. The bootstrap type `"stud"` is not supported, as the
#' statistic returns no per-replicate variance estimate.
#'
#' Note that for paired data the difference of the medians is not the
#' median of the differences. For the latter, use `medianCI(x - y)`.
#' Neither is it the Hodges-Lehmann shift estimate reported by
#' [wilcox.test()], which coincides with the median difference only under
#' a pure location shift.
#'
#' @param x a (non-empty) numeric vector of data values.
#' @param y a (non-empty) numeric vector of data values.
#' @param conf.level confidence level of the interval.
#' @param sides a character string specifying the side of the confidence
#'   interval, must be one of `"two.sided"` (default), `"left"` or
#'   `"right"`. You can specify just the initial letter. `"left"` would
#'   be analogue to a hypothesis of `"greater"` in a `t.test`.
#' @param method a character string, either `"bonett-price"` (default) or
#'   `"boot"`. The bootstrap type can be provided by the dots, see
#'   [boot::boot.ci()].
#' @param paired a logical indicating whether `x` and `y` are paired.
#'   Defaults to `FALSE`. Requires `method = "boot"`.
#' @param na.rm logical. Should missing values be removed? Defaults to
#'   `FALSE`. For paired data, incomplete pairs are removed.
#' @param \dots further arguments passed to the bootstrap, e.g. `R` and
#'   `type`.
#'
#' @return A named numeric vector with elements:
#' \describe{
#'   \item{`mediandiff`}{point estimate, the difference
#'     `median(x) - median(y)`}
#'   \item{`lci`}{lower confidence interval bound}
#'   \item{`uci`}{upper confidence interval bound}
#' }
#'
#' @references
#' Bonett, D. G. and Price, R. M. (2002). Statistical inference for a
#' linear function of medians: confidence intervals, hypothesis testing,
#' and sample size requirements. *Psychological Methods* **7**(3),
#' 370--383. \doi{10.1037//1082-989x.7.3.370}
#'
#' Price, R. M. and Bonett, D. G. (2002). Distribution-free confidence
#' intervals for difference and ratio of medians. *Journal of Statistical
#' Computation and Simulation* **72**(2), 119--124.
#' \doi{10.1080/00949650212140}
#'
#' @seealso [meanDiffCI()], [medianCI()], [boot::boot.ci()]
#'
#' @examples
#' x <- mtcars[mtcars$am == 0, "mpg"]
#' y <- mtcars[mtcars$am == 1, "mpg"]
#'
#' medianDiffCI(x, y)
#' medianDiffCI(x, y, conf.level = 0.99)
#' medianDiffCI(x, y, sides = "left")
#'
#' # bootstrap
#' set.seed(1)
#' medianDiffCI(x, y, method = "boot", type = "perc")
#' medianDiffCI(x, y, method = "boot", type = "bca")
#'
#' # paired: difference of medians vs. median of differences
#' set.seed(2)
#' drug2 <- sleep$extra[sleep$group == 2]
#' drug1 <- sleep$extra[sleep$group == 1]
#' medianDiffCI(drug2, drug1, method = "boot", paired = TRUE, type = "perc")
#' medianCI(drug2 - drug1)
#'
#' @family ci.location
#' @concept confidence-interval
#' @concept median
#'
#' @export
medianDiffCI <- function(x,
                         y,
                         conf.level = 0.95,
                         sides = c("two.sided", "left", "right"),
                         method = c("bonett-price", "boot"),
                         paired = FALSE,
                         na.rm = FALSE,
                         ...) {

  if (!is.numeric(x) || !is.numeric(y))
    stop("'x' and 'y' must be numeric")

  if (paired && length(x) != length(y))
    stop("For paired intervals, 'x' and 'y' must have equal length")

  method <- match.arg(method)
  sides  <- match.arg(sides)

  if (paired && method != "boot")
    stop("paired intervals require method = \"boot\"")

  if (na.rm) {
    if (paired) {
      ok <- complete.cases(x, y)
      x <- x[ok]
      y <- y[ok]
    } else {
      x <- x[!is.na(x)]
      y <- y[!is.na(y)]
    }
  } else if (anyNA(x) || anyNA(y)) {
    return(c(mediandiff = NA_real_, lci = NA_real_, uci = NA_real_))
  }

  if (length(x) < 2 || length(y) < 2)
    stop("Need at least two observations per group")

  .checkSidedLevel(conf.level, sides)

  if (sides != "two.sided")
    conf.level <- 1 - 2 * (1 - conf.level)

  res <- switch(
    method,
    "bonett-price" = .medianDiffCI.bp(x, y, conf.level = conf.level),
    boot           = .medianDiffCI.boot(x, y, conf.level = conf.level,
                                        paired = paired, ...)
  )

  if (sides == "left")
    res["uci"] <- Inf
  else if (sides == "right")
    res["lci"] <- -Inf

  res
}



# == internal helper functions ================================================

#' @keywords internal
.medianSE.bp <- function(x) {

  # Bonett & Price (2002): SE of a sample median from the order statistics
  # Y(c) and Y(n-c+1), scaled by the normal quantile that matches the
  # exact binomial coverage of that pair
  n <- length(x)
  x <- sort(x)
  # c as in the authors' own implementation (statpsych::ci.median2),
  # clamped to 1 for n <= 4, where round() would give 0
  c <- max(1L, as.integer(round(n / 2 - sqrt(n))))
  z <- qnorm(pbinom(c - 1L, n, 0.5), lower.tail = FALSE)

  (x[n - c + 1L] - x[c]) / (2 * z)
}


#' @keywords internal
.medianDiffCI.bp <- function(x, y, conf.level) {

  est <- median(x) - median(y)
  se  <- sqrt(.medianSE.bp(x)^2 + .medianSE.bp(y)^2)
  z   <- qnorm(1 - (1 - conf.level) / 2)

  c(mediandiff = est, lci = est - z * se, uci = est + z * se)
}


#' @keywords internal
.medianDiffCI.boot <- function(x, y, conf.level, paired = FALSE, ...) {

  args <- .extractBootArgs(list(...))

  if (args$type == "stud")
    stop("bootstrap type 'stud' is not supported by medianDiffCI(): ",
         "the statistic function does not return a per-replicate ",
         "variance estimate; use type = 'perc', 'bca', 'basic' or 'norm'")

  if (paired) {

    # resample whole pairs, keeping the dependence between x and y
    boot.fun <- boot::boot(
      cbind(x, y),
      statistic = function(d, i) median(d[i, 1L]) - median(d[i, 2L]),
      R         = args$R,
      parallel  = args$parallel,
      ncpus     = args$ncpus
    )

  } else {

    # stratified resampling: each replicate keeps both group sizes
    g <- rep(1:2, c(length(x), length(y)))

    boot.fun <- boot::boot(
      c(x, y),
      statistic = function(d, i) {
        gi <- g[i]
        median(d[i][gi == 1L]) - median(d[i][gi == 2L])
      },
      strata    = g,
      R         = args$R,
      parallel  = args$parallel,
      ncpus     = args$ncpus
    )
  }

  ci <- boot::boot.ci(boot.fun, conf = conf.level, type = args$type)

  # by name, not ci[[4]]: boot.ci() drops a component it cannot compute
  bnd <- .bootCIBounds(ci, args$type)

  c(mediandiff = boot.fun$t0, lci = bnd[1L], uci = bnd[2L])
}
