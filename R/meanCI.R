
#' Confidence Intervals for the Mean
#'
#' Collection of several approaches to determine confidence intervals for the
#' mean. Both, the classical way and bootstrap intervals are implemented for
#' both, normal and trimmed means.
#'
#' The confidence intervals for the trimmed means use winsorized variances as
#' described in the references.
#'
#' The bootstrap type \code{"stud"} (studentized) requires a variance
#' estimate to be returned alongside the point estimate on every bootstrap
#' replicate; this is supported for both the trimmed and untrimmed mean.
#'
#' @param x a (non-empty) numeric vector of data values.
#' @param conf.level confidence level of the interval.
#' @param sides a character string specifying the side of the confidence
#' interval, must be one of \code{"two.sided"} (default), \code{"left"} or
#' \code{"right"}. \code{"left"} would be analogue to a hypothesis of
#' \code{"greater"} in a \code{t.test}. You can specify just the initial
#' letter.
#' @param method A vector of character strings representing the type of
#' intervals required. The value should be any subset of the values
#' \code{"classic"}, \code{"boot"}.  See \code{\link[boot]{boot.ci}}.
#' @param sd the standard deviation of x. If provided it's interpreted as sd of
#' the population and the normal quantiles will be used for constructing the
#' confidence intervals. If left to \code{NULL} (default) the sample
#' \code{sd(x)} will be calculated and used in combination with the
#' t-distribution.
#' @param trim the fraction (0 to 0.5) of observations to be trimmed from each
#' end of \code{x} before the mean is computed. Values of \code{trim} outside
#' that range are taken as the nearest endpoint.
#' @param na.rm a logical value indicating whether \code{NA} values should be
#' stripped before the computation proceeds. Defaults to FALSE.
#'
#' @param ... further arguments are passed to the \code{\link[boot]{boot}} function.
#' Supported arguments are \code{type} (\code{"norm"}, \code{"basic"},
#' \code{"stud"}, \code{"perc"}, \code{"bca"}), \code{parallel} and the number
#' of bootstrap replicates \code{R}. If not defined those will be set to their
#' defaults, being \code{"basic"} for \code{type}, option
#' \code{"boot.parallel"} (and if that is not set, \code{"no"}) for
#' \code{parallel} and \code{999} for \code{R}.
#'
#' @return A named numeric vector with elements:
#' \describe{
#'   \item{\code{est}}{point estimate}
#'   \item{\code{lci}}{lower confidence interval bound}
#'   \item{\code{uci}}{upper confidence interval bound}
#' }
#'
#' @seealso \code{DescToolsX::meanX}, \code{\link{t.test}}, \code{\link{meanDiffCI}},
#' \code{\link{medianCI}}, \code{\link{varCI}}, \code{\link{meanCIn}}
#'
#' @references Wilcox, R. R., Keselman H. J. (2003) Modern robust data analysis
#' methods: measures of central tendency \emph{Psychol Methods}, 8(3):254-74
#'
#' Wilcox, R. R. (2005) \emph{Introduction to robust estimation and hypothesis
#' testing} Elsevier Academic Press
#'
#' @examples
#'
#' x <- mtcars$mpg
#'
#' meanCI(x, na.rm=TRUE)
#' meanCI(x, conf.level=0.99, na.rm=TRUE)
#'
#' meanCI(x, sides="left", na.rm=TRUE)
#' # same as:
#' t.test(x, alternative="greater")
#'
#' meanCI(x, sd=25, na.rm=TRUE)
#'
#' # the different types of bootstrap confints
#' meanCI(x, method="boot", type="norm", na.rm=TRUE)
#' meanCI(x, trim=0.1, method="boot", type="norm", na.rm=TRUE)
#' meanCI(x, trim=0.1, method="boot", type="basic", na.rm=TRUE)
#' meanCI(x, trim=0.1, method="boot", type="stud", na.rm=TRUE)
#' meanCI(x, trim=0.1, method="boot", type="perc", na.rm=TRUE)
#' meanCI(x, trim=0.1, method="boot", type="bca", na.rm=TRUE)
#'
#' meanCI(x, trim=0.1, method="boot", type="bca", R=1999, na.rm=TRUE)
#'
#' # Getting the meanCI for more than 1 column
#' round(t(sapply(mtcars[, c("mpg", "hp")], meanCI, na.rm=TRUE)), 3)
#'
#' @family ci.location
#' @concept confidence-interval
#'
#' @export
meanCI <- function(x,
                   conf.level = 0.95,
                   sides = c("two.sided", "left", "right"),
                   method = c("classic", "boot"),
                   sd = NULL,
                   trim = 0,
                   na.rm = FALSE,
                   ...) {

  if (na.rm)
    x <- na.omit(x)

  if (!is.numeric(x))
    stop("'x' must be numeric")

  if (length(x) < 2)
    stop("Need at least two observations")

  if (!is.numeric(trim) || length(trim) != 1 ||
      trim < 0 || trim >= 0.5)
    stop("'trim' must be a single number in [0, 0.5)")

  method <- match.arg(method)

  sides <- match.arg(
    sides,
    choices = c("two.sided", "left", "right"),
    several.ok = FALSE
  )

  if (sides != "two.sided")
    conf.level <- 1 - 2 * (1 - conf.level)

  res <- switch(

    method,

    classic = .meanCI.classic(
      x,
      conf.level = conf.level,
      sd = sd,
      trim = trim
    ),

    boot = .meanCI.boot(
      x,
      conf.level = conf.level,
      trim = trim,
      ...
    )
  )

  names(res) <- c("est", "lci", "uci")

  if (sides == "left") {

    res["uci"] <- Inf

  } else if (sides == "right") {

    res["lci"] <- -Inf
  }

  res
}



# == internal helper functions ================================================


# calculate winsorized variance and corresponding degrees of freedom
#' @keywords internal
.winvar <- function(x, trim) {

  n <- length(x)

  # calculate the winsorized variance of x
  trn <- floor(trim * n) + 1

  # new 17.2.2015:
  minval <- sort(x, partial = trn)[trn]
  maxval <- sort(x, partial = max((n - trn + 1), 1))[max((n - trn + 1), 1)]

  wvar <- var(winsorize(x, val = c(minval, maxval)))

  # This was an overkill, we need only the n-thest value here:
  # winvar <- var(winsorize(x,
  #   minval=max(Small(x, trn)),
  #   maxval=min(Large(x, trn)))
  # )

  # degrees of freedom
  DF <- n - 2 * (trn - 1) - 1

  c(var = wvar, DF = DF)
}



#' @keywords internal
.meanCI.classic <- function(x,
                            conf.level,
                            sd = NULL,
                            trim = 0) {

  if (trim != 0) {

    # see:
    # http://dornsife.usc.edu/assets/sites/239/docs/Rallfun-v27.txt
    # http://www.psychology.mcmaster.ca/bennett/boot09/rt2.pdf

    wvar <- .winvar(x, trim)

    # the standard error
    se <- sqrt(wvar["var"]) /
      ((1 - 2 * trim) * sqrt(length(x)))

    mean(x, trim = trim) +
      c(0, -1, 1) *
      qt(
        1 - (1 - conf.level) / 2,
        wvar["DF"]
      ) * se

  } else {

    if (is.null(sd)) {

      a <- qt(
        p = (1 - conf.level) / 2,
        df = length(x) - 1
      ) * sd(x) / sqrt(length(x))

    } else {

      a <- qnorm(
        p = (1 - conf.level) / 2
      ) * sd / sqrt(length(x))
    }

    c(
      mean(x),
      mean(x) + a,
      mean(x) - a
    )
  }
}



#' @keywords internal
.meanCI.boot <- function(x,
                         conf.level,
                         trim = 0,
                         ...) {

  # see:
  # http://www.psychology.mcmaster.ca/bennett/boot09/percentileT.pdf

  args <- .extractBootArgs(list(...))

  # we need separate functions for trimmed means and normal means
  if (trim != 0) {

    boot.fun <- boot::boot(

      x,

      function(x, i) {

        # this is according to the example in boot.ci, with the
        # winsorized variance recomputed on the RESAMPLED data x[i]
        # (not the original x) so it varies correctly across
        # replicates, and subset to its 'var' component only, so the
        # statistic vector returned here always has exactly 2 elements
        m <- mean(
          x[i],
          na.rm = FALSE,
          trim = trim
        )

        v <- unname(.winvar(x[i], trim)["var"]) /
          ((1 - 2 * trim) * sqrt(length(i)))^2

        c(m, v)
      },

      R        = args$R,
      parallel = args$parallel,
      ncpus    = args$ncpus
    )

  } else {

    boot.fun <- boot::boot(

      x,

      function(x, i) {

        # this is according to the example in boot.ci
        m <- mean(
          x[i],
          na.rm = FALSE
        )

        n <- length(i)

        # variance of the sample mean: var(x[i]) / n
        v <- var(x[i]) / n

        c(m, v)

        # IMPORTANT:
        # boot.ci requires the estimated VARIANCE
        # of the statistic
      },

      R        = args$R,
      parallel = args$parallel,
      ncpus    = args$ncpus
    )
  }

  ci <- boot::boot.ci(
    boot.fun,
    conf = conf.level,
    type = args$type
  )

  if (args$type == "norm") {

    c(
      mean = boot.fun$t0[1],
      lci  = ci[[4]][2],
      uci  = ci[[4]][3]
    )

  } else {

    c(
      mean = boot.fun$t0[1],
      lci  = ci[[4]][4],
      uci  = ci[[4]][5]
    )
  }
}
