
#' Confidence Interval for a Variance
#'
#' Computes confidence intervals for a population variance using
#' classical chi-square, Bonett, or bootstrap methods.
#' 
#' The confidence interval for the variance is very sensitive to non-normality
#' in the data. Bonett (2006) has proposed an interval that is nearly exact
#' when the data is normally distributed and provides good performance for
#' moderately non-normal data. See the references for the details.
#' 
#' @param x a (non-empty) numeric vector of data values.
#' @param conf.level confidence level of the interval.
#' @param sides a character string specifying the side of the confidence
#' interval, must be one of \code{"two.sided"} (default), \code{"left"} or
#' \code{"right"}.  You can specify just the initial letter. \code{"left"}
#' would be analogue to a hypothesis of \code{"greater"} in a \code{t.test}.
#' @param method vector of character strings representing the type of intervals
#' required.  The value should be any subset of the values \code{"classic"},
#' \code{"bonett"}, \code{"norm"}, \code{"boot"}. Bootstrap type can be set by
#' the ... arguments.  See \code{\link[boot]{boot.ci}}.
#' @param na.rm logical. Should missing values be removed? Defaults to FALSE.
#' @param \dots further arguments, can be used to provide further arguments to
#' the boot function.
#' 
#' @return Named numeric vector with:
#' \itemize{
#'   \item \code{var}: sample variance
#'   \item \code{lci}: lower confidence limit
#'   \item \code{uci}: upper confidence limit
#' }
#' @seealso \code{\link{meanCI}}, \code{\link{medianCI}},
#' \code{\link{varTest}}, \code{\link[DescToolsX]{varX}} 
#' 
#' @references Bonett (2006) Approximate Confidence Interval for Standard
#' Deviation of Nonnormal Distributions, \emph{Computational Statistics and
#' Data Analysis}, Vol. 50, pp. 775 - 782.\cr
#' https://www.itl.nist.gov/div898/software/dataplot/refman1/auxillar/sdconfli.htm
#' (might be outdated)
#' @keywords univar
#' @examples
#' 
#' x <- mtcars$mpg
#' 
#' varCI(x, na.rm=TRUE)
#' varCI(x, conf.level=0.99, na.rm=TRUE)
#' 
#' x <- c(14.816, 14.863, 14.814, 14.998, 14.965, 14.824, 14.884, 14.838, 14.916,
#'        15.021, 14.874, 14.856, 14.860, 14.772, 14.980, 14.919)
#' varCI(x, conf.level=0.9)
#' 
#' # and for the standard deviation
#' sqrt(varCI(x, conf.level=0.9))
#' 
#' 
#' # from Bonett's paper
#' # expected results:
#' # ------------------------------------
#' #  conf.lvl       sd      lci      uci
#' # ------------------------------------
#' #      90.0   0.5168   0.3592   0.9359
#' #      95.0   0.5168   0.3263   1.0841
#' #      99.0   0.5168   0.2607   1.5109
#' 
#' p <- c(15.83, 16.01, 16.24, 16.42, 15.33, 15.44, 16.88, 16.31)
#' sqrt(varCI(p, method="bonett", conf.level=0.9))
#' sqrt(varCI(p, method="bonett"))
#' sqrt(varCI(p, method="bonett", conf.level=0.99))
#' 
#' # some bootstrap intervals
#' varCI(x, method="boot", type="norm")
#' varCI(x, method="boot", type="perc")
#' varCI(x, method="boot", type="bca")
#' 




#' @family test.variance  
#' @concept variance-test  
#' @concept confidence-interval
#'
#'
#' @export
varCI <- function(x,
                  conf.level = 0.95,
                  sides = c("two.sided", "left", "right"),
                  method = c("classic", "bonett", "boot"),
                  na.rm = FALSE,
                  ...) {
  
  if (na.rm)
    x <- na.omit(x)
  
  if (!is.numeric(x))
    stop("'x' must be numeric")
  
  if (length(x) < 2)
    stop("Need at least two observations")
  
  method <- match.arg(method)
  sides  <- match.arg(sides)
  
  if (sides != "two.sided")
    conf.level <- 1 - 2 * (1 - conf.level)
  
  res <- switch(method,
    classic = .varCI.classic(x, conf.level = conf.level),
    bonett  = .varCI.bonett(x, conf.level = conf.level),
    boot    = .varCI.boot(x, conf.level = conf.level, ...)
  )
  
  if (sides == "left") {
    res["uci"] <- Inf
  } else if (sides == "right") {
    res["lci"] <- 0
  }
  
  res
}


# == internal helper functions ===============================================

.varCI.classic <- function(x, conf.level) {
  
  df <- length(x) - 1
  v  <- var(x)
  
  c(var = v,
    lci = df * v / qchisq( (1 - conf.level) / 2, df, lower.tail = FALSE),
    uci = df * v / qchisq((1 - conf.level) / 2, df)
  )
}


.varCI.bonett <- function(x, conf.level) {
  
  n <- length(x)
  
  if (n <= 4)
    stop("Bonett method requires n > 4")
  
  z <- qnorm(1 - (1 - conf.level) / 2)
  
  cc <- n / (n - z)
  
  v   <- var(x)
  mtr <- mean(x, trim = 1 / (2 * sqrt(n - 4)))
  m   <- mean(x)
  
  gam4 <- n * sum((x - mtr)^4) / (sum((x - m)^2))^2
  
  se <- cc * sqrt((gam4 - (n - 3) / n) / (n - 1))
  
  lci <- exp(log(cc * v) - z * se)
  uci <- exp(log(cc * v) + z * se)
  
  c(var = v,
    lci = lci,
    uci = uci
  )
}



.varCI.boot <- function(x, conf.level, ...) {
  
  args <- .extractBootArgs(list(...))
  
  boot.fun <- boot::boot(
    x,
    statistic = function(x, d)
      var(x[d]),
    R        = args$R,
    parallel = args$parallel,
    ncpus    = args$ncpus
  )
  
  ci <- boot::boot.ci(
    boot.fun,
    conf = conf.level,
    type = args$type
  )
  
  if (args$type == "norm") {
    
    c(
      var = boot.fun$t0,
      lci = ci[[4]][2],
      uci = ci[[4]][3]
    )
    
  } else {
    
    c(
      var = boot.fun$t0,
      lci = ci[[4]][4],
      uci = ci[[4]][5]
    )
  }
}
