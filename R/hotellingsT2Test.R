
#' Hotelling's T2 Test
#' 
#' A multivariate generalization of the two-sample t-test, testing whether 
#' two groups differ significantly across multiple dependent variables 
#' simultaneously.
#' 
#' 
#' Hotelling's T2 test is the multivariate generlisation of the Student's t
#' test. A one-sample Hotelling's T2 test can be used to test if a set of
#' vectors of data (which should be a sample of a single statistical
#' population) has a mean equal to a hypothetical mean. A two-sample
#' Hotelling's T2 test may be used to test for significant differences between
#' the mean vectors (multivariate means) of two multivariate data sets are
#' different.
#' 
#' 
#' The classical test for testing the location of a multivariate population or
#' for testing the mean difference for two multivariate populations. When
#' \code{test = "f"} the F-distribution is used for the test statistic and it
#' is assumed that the data are normally distributed. If the chisquare
#' approximation is used, the normal assumption can be relaxed to existence of
#' second moments.  In the two sample case both populations are assumed to have
#' the same covariance matrix.
#' 
#' The formula interface is only applicable for the 2-sample tests.
#' 
#' @name hotellingsT2Test
#' @aliases hotellingsT2Test hotellingsT2Test.default hotellingsT2Test.formula
#' @param x a numeric data frame or matrix.
#' @param y an optional numeric data frame or matrix for the two sample test.
#' If \code{NULL} a one sample test is performed.
#' @param mu a vector indicating the hypothesized value of the mean (or
#' difference in means if a two sample test is performed). \code{NULL}
#' represents origin or no difference between the groups.
#' @param test if \code{"f"}, the decision is based on the F-distribution, if
#' \code{"chi"} a chi-squared approximation is used.
#' @param formula a formula of the form \code{x ~ g} where \code{x} is a
#' numeric matrix giving the data values and \code{g} a factor with two levels
#' giving the corresponding groups.
#' @param data an optional matrix or data frame (or similar: see
#' \code{\link{model.frame}}) containing the variables in the formula
#' \code{formula}.  By default the variables are taken from
#' \code{environment(formula)}.
#' @param subset an optional vector specifying a subset of observations to be
#' used.
#' @param na.action a function which indicates what should happen when the data
#' contain NAs. Defaults to \code{getOption("na.action")}.
#' @param \dots further arguments to be passed to or from methods.
#' @return A list with class 'htest' containing the following components:
#' \item{statistic }{the value of the T2-statistic. (That is the scaled value
#' of the statistic that has an F distribution or a chisquare distribution
#' depending on the value of \code{test}).} \item{parameter}{the degrees of
#' freedom for the T2-statistic.} \item{p.value}{the p-value for the test.}
#' \item{null.value}{the specified hypothesized value of the mean or mean
#' difference depending on whether it was a one-sample test or a two-sample
#' test.} \item{alternative}{a character string with the value 'two.sided'.}
#' \item{method}{a character string indicating what type of test was
#' performed.} \item{data.name}{a character string giving the name of the data
#' (and grouping vector).}
#' 
#' @note
#' Adapted from code by Klaus Nordhausen to conform to package standards.
#' 
#' @references Nordhausen K., Sirkia S., Oja H. and Tyler D. E. (2012)
#' \emph{ICSNP: Tools for Multivariate Nonparametrics}. R package version
#' 1.0-9.\cr \url{https://cran.r-project.org/package=ICSNP}
#' 
#' Anderson, T.W. (2003), \emph{An introduction to multivariate analysis}, New
#' Jersey: Wiley.
#' 
#' @examples
#' 
#' math.teach <- data.frame(
#'   teacher = factor(rep(1:2, c(3, 6))),
#'   satis   = c(1, 3, 2, 4, 6, 6, 5, 5, 4),
#'   know    = c(3, 7, 2, 6, 8, 8, 10, 10, 6))
#' 
#' with(math.teach,
#'   hotellingsT2Test(cbind(satis, know) ~ teacher))
#' 

#' @rdname hotellingsT2Test
#' @family test.location
#' @concept hypothesis-testing
#' @concept multivariate
#'
#'
#' @export
hotellingsT2Test <- function(x,...) {
  UseMethod("hotellingsT2Test")
}

#' @rdname hotellingsT2Test
#' @export
hotellingsT2Test.default <- function(x, y=NULL, mu=NULL, test="f",...) {
  
  
  `HotellingsT.internal`  <-  function(x, y=NULL, mu, test) {
    n <- dim(x)[1]
    p <- dim(x)[2]
    
    if(is.null(y))     #one sample case
    {
      test.statistic <- n*as.numeric(t(colMeans(x)-mu)%*%solve(cov(x))%*%(colMeans(x)-mu))*switch(test,f=(n-p)/(p*(n-1)),chi=1)
      df.1 <- p
      df.2 <- switch(test,f=n-p,chi=NA)
      p.value <- 1-switch(test,f=pf(test.statistic,df.1,df.2),chi=pchisq(test.statistic,df.1))
      return(list(test.statistic=test.statistic,p.value=p.value,df.1=df.1,df.2=df.2))
    }
    
    # else two sample case
    n1 <- n
    n2 <- dim(y)[1]
    xmeans <- colMeans(x)
    ymeans <- colMeans(y)
    x.diff <- sweep(x,2,xmeans)
    y.diff <- sweep(y,2,ymeans)
    S.pooled <- 1/(n1+n2-2)*(t(x.diff)%*%x.diff+t(y.diff)%*%y.diff)
    test.statistic <- n1*n2/(n1+n2)*t(xmeans-ymeans-mu)%*%solve(S.pooled)%*%(xmeans-ymeans-mu)*switch(test,f=(n1+n2-p-1)/(p*(n1+n2-2)),chi=1)
    df.1 <- p
    df.2 <- switch(test,f=n1+n2-p-1,chi=NA)
    p.value <- 1-switch(test,f=pf(test.statistic,df.1,df.2),chi=pchisq(test.statistic,df.1))
    list(test.statistic=test.statistic,p.value=p.value,df.1=df.1,df.2=df.2)
  }
  
  
  if (is.null(y)) {
    DNAME <- deparse(substitute(x))
  } else {
    DNAME=paste(deparse(substitute(x)),"and",deparse(substitute(y)))
  }
  
  xok <- complete.cases(x)
  x <- x[xok,]
  if(!all(sapply(x, is.numeric))) stop("'x' must be numeric")
  x <- as.matrix(x)
  
  p <- dim(x)[2]
  
  if (!is.null(y)) {
    yok <- complete.cases(y)
    y <- y[yok,]
    
    if(!all(sapply(y, is.numeric))) stop("'y' must be numeric")
    if (p!=dim(y)[2]) stop("'x' and 'y' must have the same number of columns")
    y <- as.matrix(y)
  }
  
  if (is.null(mu)) mu <- rep(0,p)
  else if (length(mu)!=p) stop("length of 'mu' must equal the number of columns of 'x'")
  
  test <- match.arg(test,c("f","chi"))
  
  if (is.null(y) & test=="f") version <- "one.sample.f"
  if (is.null(y) & test=="chi") version <- "one.sample.chi"
  if (!is.null(y) & test=="f") version <- "two.sample.f"
  if (!is.null(y) & test=="chi") version <- "two.sample.chi"
  
  res1 <- switch(version,
                 "one.sample.f"={
                   result <- HotellingsT.internal(x,mu=mu,test=test)
                   STATISTIC <- result$test.statistic
                   names(STATISTIC) <- "T.2"
                   PVAL <- result$p.value
                   METHOD <- "Hotelling's one sample T2-test"
                   PARAMETER <- c(result$df.1,result$df.2)
                   names(PARAMETER) <- c("df1","df2")
                   RVAL <- list(statistic=STATISTIC,p.value=PVAL,method=METHOD,parameter=PARAMETER)
                   
                   RVAL}
                 ,
                 "one.sample.chi"={
                   result <- HotellingsT.internal(x,mu=mu,test=test)
                   STATISTIC <- result$test.statistic
                   names(STATISTIC) <- "T.2"
                   PVAL <- result$p.value
                   METHOD <- "Hotelling's one sample T2-test"
                   PARAMETER <- c(result$df.1)
                   names(PARAMETER) <- c("df")
                   RVAL <- list(statistic=STATISTIC,p.value=PVAL,method=METHOD,parameter=PARAMETER)
                   
                   RVAL}
                 ,
                 "two.sample.f"={
                   result <- HotellingsT.internal(x,y,mu,test)
                   STATISTIC <- result$test.statistic
                   names(STATISTIC) <- "T.2"
                   PVAL <- result$p.value
                   METHOD <- "Hotelling's two sample T2-test"
                   PARAMETER <- c(result$df.1,result$df.2)
                   names(PARAMETER) <- c("df1","df2")
                   RVAL <- list(statistic=STATISTIC,p.value=PVAL,method=METHOD,parameter=PARAMETER)
                   
                   RVAL}
                 ,
                 "two.sample.chi"={
                   result <- HotellingsT.internal(x,y,mu,test)
                   STATISTIC <- result$test.statistic
                   names(STATISTIC) <- "T.2"
                   PVAL <- result$p.value
                   METHOD <- "Hotelling's two sample T2-test"
                   PARAMETER <- c(result$df.1)
                   names(PARAMETER) <- c("df")
                   RVAL <- list(statistic=STATISTIC,p.value=PVAL,method=METHOD,parameter=PARAMETER)
                   
                   RVAL}
  )
  ALTERNATIVE="two.sided"
  NVAL <- paste("c(",paste(mu,collapse=","),")",sep="")
  if (is.null(y)) names(NVAL) <- "location" else names(NVAL) <- "location difference"
  res <- c(res1,list(data.name=DNAME,alternative=ALTERNATIVE,null.value=NVAL))
  class(res) <- "htest"
  return(res)
}



#' @rdname hotellingsT2Test
#' @export
hotellingsT2Test.formula <- function (formula, data, subset, na.action, ...) {
  
  if (missing(formula) || (length(formula) != 3L) || (length(attr(terms(formula[-2L]),
                                                                  "term.labels")) != 1L))
    stop("'formula' missing or incorrect")
  m <- match.call(expand.dots = FALSE)
  if (is.matrix(eval(m$data, parent.frame())))
    m$data <- as.data.frame(data)
  m[[1L]] <- quote(stats::model.frame)
  m$... <- NULL
  mf <- eval(m, parent.frame())
  DNAME <- paste(names(mf), collapse = " by ")
  names(mf) <- NULL
  response <- attr(attr(mf, "terms"), "response")
  g <- factor(mf[[-response]])
  if (nlevels(g) != 2L)
    stop("grouping factor must have exactly 2 levels")
  # DATA <- setNames(split(mf[[response]], g), c("x", "y"))
  DATA <- setNames(split(as.data.frame(mf[[response]]), g), c("x", "y"))
  
  y <- do.call("hotellingsT2Test", c(DATA, list(...)))
  y$data.name <- DNAME
  y
}


