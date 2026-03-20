
# Confidence intervall for the median

# this is used in signTest and would require a recursive dependency
# thus it's copied in here as internal function


.medianCI <- function(x, 
                     conf.level=0.95, sides = c("two.sided","left","right"), 
                     method=c("exact","boot"),
                     na.rm=FALSE, ...) {
  
  if(na.rm) x <- na.omit(x)
  
  medianCI_Binom <- function( x, conf.level = 0.95,
                              sides = c("two.sided", "left", "right"), na.rm = FALSE ){
    
    # http://www.stat.umn.edu/geyer/old03/5102/notes/rank.pdf
    # http://de.scribd.com/doc/75941305/Confidence-Interval-for-Median-Based-on-Sign-Test
    if(na.rm) x <- na.omit(x)
    n <- length(x)
    switch( match.arg(sides)
            , "two.sided" = {
              k <- qbinom(p = (1 - conf.level) / 2, size=n, prob=0.5, lower.tail=TRUE)
              ci <- sort(x)[c(k, n - k + 1)]
              attr(ci, "conf.level") <- 1 - 2 * pbinom(k-1, size=n, prob=0.5)
            }
            , "left" = {
              k <- qbinom(p = (1 - conf.level), size=n, prob=0.5, lower.tail=TRUE)
              ci <- c(sort(x)[k], Inf)
              attr(ci, "conf.level") <- 1 - pbinom(k-1, size=n, prob=0.5)
            }
            , "right" = {
              k <- qbinom(p = conf.level, size=n, prob=0.5, lower.tail=TRUE)
              ci <- c(-Inf, sort(x)[k])
              attr(ci, "conf.level") <- pbinom(k, size=n, prob=0.5)
            }
    )
    # confints for small samples can be outside the observed range e.g. n < 6
    if(identical(.stripAttr(ci), NA_real_)) {
      ci <- c(-Inf, Inf)
      attr(ci, "conf.level") <- 1
    }  
    return(ci)
  }
  
  medianCI_Boot <- function(x, conf.level=0.95, sides = c("two.sided", "left", "right"), 
                            na.rm=FALSE, ...){
    
    if(sides!="two.sided")
      conf.level <- 1 - 2*(1-conf.level)
    
    R <- .inDots(..., arg="R", default=999)
    boot.med <- boot::boot(x, function(x, d) {
      median(x[d], na.rm=na.rm)
      # standard error for the median required for studentized bci type:
      # not implemented here, as not suitable for this case.
      # sqrt(pi/2) * MeanSE(x[d])
      # mad(x[d], na.rm=na.rm) / sqrt(length(na.omit(x[d])))
      
    }, R=R)
    
    dots <- list(...)
    if(is.null(dots[["type"]]))
      dots$type <- "perc"
    
    if(!(dots$type %in% c("norm","basic","perc","bca"))){
      warning(gettextf("bootstrap type '%s' is not supported", dots$type))
      return( c(NA, NA))
    }
    
    dots$boot.out <- boot.med
    dots$conf <- conf.level
    
    res <- do.call(boot::boot.ci, dots)
    
    if(dots$type == "norm")
      # uses different structure for results
      res <- res[[4]][c(2,3)]
    else
      res <- res[[4]][c(4,5)]
    
    return(res)
  }
  
  
  
  sides <- match.arg(sides, choices = c("two.sided","left","right"), several.ok = FALSE)
  
  # if(sides!="two.sided")
  #   conf.level <- 1 - 2*(1-conf.level)
  
  # alte Version, ziemlich grosse Unterschiede zu wilcox.test:
  # Bosch: Formelsammlung Statistik (bei Markus Naepflin), S. 95
  # x <- sort(x)
  # return( c(
  # x[ qbinom(alpha/2,length(x),0.5) ], ### lower limit
  # x[ qbinom(1-alpha/2,length(x),0.5) ] ### upper limit
  # ) )
  
  method <- match.arg(arg=method, choices=c("exact","boot"))
  
  switch( method
          , "exact" = { # this is the SAS-way to do it
            # https://stat.ethz.ch/pipermail/r-help/2003-September/039636.html
            r <- medianCI_Binom(x, conf.level = conf.level, sides=sides)
          }
          , "boot" = {
            r <- medianCI_Boot(x, conf.level = conf.level, sides=sides, ...)
          } )
  
  med <- median(x, na.rm=na.rm)
  if(is.na(med)) {   # do not report a CI if the median is not defined...
    res <- rep(NA, 3)
    
  } else {
    res <- c(median=med, r)
    # report the conf.level which can deviate from the required one
    if(method=="exact")  attr(res, "conf.level") <-  attr(r, "conf.level")
  }
  names(res) <- c("median","lwr.ci","upr.ci")
  
  if(sides=="left")
    res[3] <- Inf
  else if(sides=="right")
    res[2] <- -Inf
  
  return( res )
  
}


# == internal helper functions ========================================


.inDots <- function(..., arg, default){
  
  # was arg in the dots-args? parse dots.arguments
  arg <- unlist(match.call(expand.dots=FALSE)$...[arg])
  
  # if arg was not in ... then return default
  if(is.null(arg)) arg <- default
  
  return(arg)
  
}



.stripAttr <- function(x, attr_names=NULL) {
  
  if(is.null(attr_names))
    attributes(x) <- NULL
  else
    for(a in attr_names) 
      attr(x, which = a) <- NULL
    
    return(x)
}





