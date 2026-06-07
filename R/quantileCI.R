
#' Confidence Interval for Any Quantile 
#' 
#' Calculates the confidence interval for any quantile. Although bootstrapping
#' might be a good approach for getting senisble confidence intervals there's
#' sometimes need to have a nonparameteric alternative. This function offers
#' one. 
#' 
#' The \code{"exact"} method corresponds to the way the confidence interval for
#' the median is calculated in SAS. \cr The boot confidence interval type is
#' calculated by means of \code{\link[boot]{boot.ci}} with default type
#' \code{"basic"}. 
#' 
#' @param x a (non-empty) numeric vector of data values.
#' @param conf.level confidence level of the interval
#' @param sides a character string specifying the side of the confidence
#' interval, must be one of \code{"two.sided"} (default), \code{"left"} or
#' \code{"right"} (abbreviations allowed). \cr\code{"left"} would be analogue
#' to a \code{"greater"} hypothesis in a \code{t.test}.
#' @param method defining the type of interval that should be calculated (one
#' out of \code{"exact"}, \code{"boot"}). Default is \code{"exact"}. See
#' Details.
#' @param probs numeric vector of probabilities with values in \emph{\verb{[0,1]}}.
#' (Values up to \code{2e-14} outside that range are accepted and moved to the
#' nearby endpoint.)
#' @param na.rm logical. Should missing values be removed? Defaults to
#' \code{FALSE}.
#' @param \dots bootstrap arguments can be provided by the dots argument. See
#' \code{\link[boot]{boot.ci}} for details.
#' 
#' @return
#' A matrix with columns \code{est} (estimated quantile), \code{lci} (lower
#' confidence limit), and \code{uci} (upper confidence limit), with one row
#' per element of \code{probs}. For the \code{"exact"} method, an attribute
#' \code{conf.level} reports the achieved coverage (which may differ from the
#' requested level).
#' @note based on code of W Huber on StackExchange
#' @seealso \code{\link[DescToolsX]{quantileX}}, \code{\link{quantile}},
#' \code{\link{medianCI}}
#' @examples
#' 
#' x <- mtcars$mpg
#' quantileCI(x, probs=0.25, na.rm=TRUE)
#' 
#' quantileCI(x, na.rm=TRUE)
#' quantileCI(x, conf.level=0.99, na.rm=TRUE)
#' 
#' # multiple probs
#' quantileCI(1:100, method="exact" , probs = c(0.25, 0.75, .80, 0.95))
#' quantileCI(1:100, method="boot" , probs = c(0.25, 0.75, .80, 0.95))
#'
#'


#' @family ci.location
#' @concept confidence-intervals
#' @concept descriptive-statistics
#' @concept nonparametric
#'
#'

#' @export
quantileCI <- function(x,
                       conf.level = 0.95,
                       sides = c("two.sided", "left", "right"),
                       method = c("exact", "boot"),
                       probs = seq(0, 1, .25),
                       na.rm = FALSE,
                       ...) {
  
  if (na.rm)
    x <- na.omit(x)
  
  if (anyNA(x))
    stop(
      "missing values and NaN's not allowed if 'na.rm' is FALSE"
    )
  
  if (!is.numeric(x))
    stop("'x' must be numeric")
  
  if (length(x) < 2)
    stop("Need at least two observations")
  
  if (any(!is.finite(probs)) ||
      any(probs < 0 | probs > 1))
    stop("'probs' must contain values in [0, 1]")
  
  sides <- match.arg(
    sides,
    choices = c("two.sided", "left", "right"),
    several.ok = FALSE
  )
  
  method <- match.arg(
    arg = method,
    choices = c("exact", "boot")
  )
  
  r <- switch(
    
    method,
    
    exact = {
      
      rr <- lapply(
        
        probs,
        
        function(p)
          .quantileCI.exact(
            x,
            prob = p,
            conf.level = conf.level,
            sides = sides
          )
      )
      
      coverage <- sapply(
        rr,
        function(z)
          attr(z, "conf.level")
      )
      
      rr <- do.call(rbind, rr)
      
      attr(rr, "conf.level") <- coverage
      
      rr
    },
    
    boot = {
      
      .quantileCI.boot(
        x,
        probs = probs,
        conf.level = conf.level,
        sides = sides,
        ...
      )
    }
  )
  
  qq <- quantile(
    x,
    probs = probs,
    na.rm = FALSE
  )
  
  res <- cbind(est = qq, r)
  
  colnames(res) <- c("est", "lci", "uci")
  
  # report the conf.level which can deviate from the required one
  if (method == "exact")
    attr(res, "conf.level") <- attr(r, "conf.level")
  
  res
}



# == internal helper functions ==============================================

.quantileCI.exact <- function(x,
                              prob,
                              conf.level = 0.95,
                              sides = c("two.sided", "left", "right")) {
  
  # Near-symmetric distribution-free confidence interval
  # for a quantile `q`.
  
  # https://stats.stackexchange.com/questions/99829/
  # how-to-obtain-a-confidence-interval-for-a-percentile
  
  # Search over a small range of upper and lower order
  # statistics for the closest coverage to 1-alpha
  # (but not less than it, if possible).
  
  n <- length(x)
  
  alpha <- 1 - conf.level
  
  if (sides == "two.sided") {
    
    u <- qbinom(
      p = 1 - alpha / 2,
      size = n,
      prob = prob
    ) + (-2:2) + 1
    
    l <- qbinom(
      p = alpha / 2,
      size = n,
      prob = prob
    ) + (-2:2)
    
    u[u > n] <- Inf
    l[l < 0] <- -Inf
    
    coverage <- outer(
      l,
      u,
      function(a, b)
        pbinom(b - 1, n, prob = prob) -
        pbinom(a - 1, n, prob = prob)
    )
    
    if (max(coverage) < 1 - alpha) {
      
      i <- which(coverage == max(coverage))
      
    } else {
      
      i <- which(
        coverage ==
          min(coverage[coverage >= 1 - alpha])
      )
    }
    
    # minimal difference
    i <- i[1]
    
    # order statistics and the actual coverage
    u <- rep(u, each = 5)[i]
    l <- rep(l, 5)[i]
    
    coverage <- coverage[i]
    
  } else if (sides == "left") {
    
    l <- qbinom(
      p = alpha,
      size = n,
      prob = prob
    )
    
    u <- Inf
    
    coverage <- 1 - pbinom(
      q = l - 1,
      size = n,
      prob = prob
    )
    
  } else if (sides == "right") {
    
    l <- -Inf
    
    u <- qbinom(
      p = 1 - alpha,
      size = n,
      prob = prob
    )
    
    coverage <- pbinom(
      q = u,
      size = n,
      prob = prob
    )
  }
  
  # get the values
  if (prob %notin% c(0, 1)) {
    
    s <- sort(
      x,
      partial = c(u, l)[is.finite(c(u, l))]
    )
    
  } else {
    
    s <- sort(x)
  }
  
  res <- c(
    lci = s[l],
    uci = s[u]
  )
  
  attr(res, "conf.level") <- coverage
  
  if (sides == "left") {
    
    res[2] <- Inf
    
  } else if (sides == "right") {
    
    res[1] <- -Inf
  }
  
  res
}


.quantileCI.boot <- function(x,
                             probs,
                             conf.level = 0.95,
                             sides = c("two.sided", "left", "right"),
                             ...) {
  
  if (sides != "two.sided")
    conf.level <- 1 - 2 * (1 - conf.level)
  
  args <- .extractBootArgs(list(...))
  
  t(sapply(
    
    probs,
    
    function(p) {
      
      boot.fun <- boot::boot(
        
        x,
        
        function(x, d)
          quantile(
            x[d],
            probs = p,
            na.rm = FALSE
          ),
        
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
          lci = ci[[4]][2],
          uci = ci[[4]][3]
        )
        
      } else {
        
        c(
          lci = ci[[4]][4],
          uci = ci[[4]][5]
        )
      }
    }
  ))
}
