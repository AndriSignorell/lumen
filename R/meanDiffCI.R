
#' Confidence Interval For Difference of Means
#' 
#' Calculates the confidence interval for the difference of two means either
#' the classical way or with the bootstrap approach.
#' 
#' This function collects code from two sources. The classical confidence
#' interval is calculated by means of \code{\link{t.test}}. The bootstrap
#' intervals are strongly based on the example in \code{\link[boot]{boot}}.
#' 
#' @param x a (non-empty) numeric vector of data values.
#' @param y a (non-empty) numeric vector of data values.
#' @param conf.level confidence level of the interval.
#' @param sides a character string specifying the side of the confidence
#' interval, must be one of \code{"two.sided"} (default), \code{"left"} or
#' \code{"right"}. You can specify just the initial letter. \code{"left"} would
#' be analogue to a hypothesis of \code{"greater"} in a \code{t.test}.
#' @param method a vector of character strings representing the type of
#' intervals required. The value should be any subset of the values
#' \code{"classic"}, \code{"boot"}. Bootstrap type can be provided by the dots.
#' See \code{\link[boot]{boot.ci}}.
#' @param paired a logical indicating whether you want confidence intervals for
#' a paired design. Defaults to \code{FALSE}.
#' @param var.equal a logical variable indicating whether to treat the two
#' variances as being equal. Default is \code{FALSE}. If \code{TRUE} then the
#' pooled variance is used to estimate the variance otherwise the Welch (or
#' Satterthwaite) approximation to the degrees of freedom is used. Passed on to
#' \code{\link{t.test}()}.
#' @param na.rm logical. Should missing values be removed? Defaults to
#' \code{FALSE}.
#' @param \dots further arguments, can be used to provide further arguments to
#' the boot function.
#' @return a numeric vector with 3 elements: \item{meandiff}{the difference:
#' mean(x) - mean(y)} \item{lci}{lower bound of the confidence interval}
#' \item{uci}{upper bound of the confidence interval}
#' @seealso \code{\link{meanCI}}, \code{\link{varCI}}, \code{\link{medianCI}},
#' \code{\link[boot]{boot.ci}}
#' @keywords univar
#' @examples
#' 
#' x <- mtcars[mtcars$am == 0, "mpg"]
#' y <- mtcars[mtcars$am == 1, "mpg"]
#' 
#' meanDiffCI(x, y, na.rm=TRUE)
#' meanDiffCI(x, y, conf.level=0.99, na.rm=TRUE)
#' 
#' # the different types of bootstrap confints
#' meanDiffCI(x, y, method="boot", type="norm", na.rm=TRUE)
#' meanDiffCI(x, y, method="boot", type="basic", na.rm=TRUE)
#' # meanDiffCI(x, y, method="boot", type="stud", na.rm=TRUE)
#' meanDiffCI(x, y, method="boot", type="perc", na.rm=TRUE)
#' meanDiffCI(x, y, method="boot", type="bca", na.rm=TRUE)
#' 
#' # for long form variables
#' with(mtcars, with(split(mpg, am), 
#'   meanDiffCI(`0`, `1`) )
#' )
#' 

              
#' @family ci.location
#' @concept confidence-intervals
#' @concept descriptive-statistics
#' @concept hypothesis-testing
#'
#'
#' @export
meanDiffCI <- function(x,
                       y,
                       conf.level = 0.95,
                       sides = c("two.sided", "left", "right"),
                       method = c("classic", "boot"),
                       paired = FALSE,
                       var.equal = FALSE,
                       na.rm = FALSE,
                       ...) {
  
  if (na.rm) {
    x <- na.omit(x)
    y <- na.omit(y)
  }
  
  if (!is.numeric(x) || !is.numeric(y))
    stop("'x' and 'y' must be numeric")
  
  if (paired && length(x) != length(y))
    stop("For paired intervals, 'x' and 'y' must have equal length")
  
  if (length(x) < 2 || length(y) < 2)
    stop("Need at least two observations per group")
  
  method <- match.arg(method)
  sides  <- match.arg(sides)
  
  if (sides != "two.sided")
    conf.level <- 1 - 2 * (1 - conf.level)
  
  res <- switch(
    
    method,
    
    classic = .meanDiffCI.classic(
      x,
      y,
      conf.level = conf.level,
      paired = paired,
      var.equal = var.equal
    ),
    
    boot = .meanDiffCI.boot(
      x,
      y,
      conf.level = conf.level,
      paired = paired,
      ...
    )
  )
  
  if (sides == "left") {
    
    res["uci"] <- Inf
    
  } else if (sides == "right") {
    
    res["lci"] <- -Inf
  }
  
  res
}



# == internal helper functions ================================================

.meanDiffCI.classic <- function(x,
                                y,
                                conf.level,
                                paired = FALSE,
                                var.equal = FALSE) {
  
  tt <- t.test(
    x,
    y,
    conf.level = conf.level,
    paired = paired,
    var.equal = var.equal
  )
  
  est <- if (paired) {
    mean(x - y)
  } else {
    mean(x) - mean(y)
  }
  
  c(
    meandiff = est,
    lci = tt$conf.int[1],
    uci = tt$conf.int[2]
  )
}



.meanDiffCI.boot <- function(x,
                             y,
                             conf.level,
                             paired = FALSE,
                             ...) {
  
  args <- .extractBootArgs(list(...))
  
  if (paired) {
    
    d <- x - y
    
    boot.fun <- boot::boot(
      d,
      statistic = function(d, i)
        mean(d[i]),
      stype    = "i",
      R        = args$R,
      parallel = args$parallel,
      ncpus    = args$ncpus
    )
    
  } else {
    
    diff.means <- function(d, f) {
      
      gp1 <- 1:table(as.numeric(d[, 2]))[1]
      
      m1 <- sum(d[gp1, 1] * f[gp1]) / sum(f[gp1])
      m2 <- sum(d[-gp1, 1] * f[-gp1]) / sum(f[-gp1])
      
      m1 - m2
    }
    
    m <- cbind(
      c(x, y),
      c(rep(1, length(x)), rep(2, length(y)))
    )
    
    boot.fun <- boot::boot(
      m,
      statistic = diff.means,
      stype     = "f",
      strata    = m[, 2],
      R         = args$R,
      parallel  = args$parallel,
      ncpus     = args$ncpus
    )
  }
  
  ci <- boot::boot.ci(
    boot.fun,
    conf = conf.level,
    type = args$type
  )
  
  if (args$type == "norm") {
    
    c(
      meandiff = boot.fun$t0,
      lci = ci[[4]][2],
      uci = ci[[4]][3]
    )
    
  } else {
    
    c(
      meandiff = boot.fun$t0,
      lci = ci[[4]][4],
      uci = ci[[4]][5]
    )
  }
}
