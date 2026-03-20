
#' Yuen t-Test For Trimmed Means
#' 
#' Performs one and two sample Yuen t-tests for trimmed means on vectors of
#' data.
#' 
#' 
#' @name yuenTTest
#' @aliases yuenTTest yuenTTest.default yuenTTest.formula
#' @param x numeric vector of data values. Non-finite (e.g. infinite or
#' missing) values will be omitted.
#' @param y an optional numeric vector of data values: as with x non-finite
#' values will be omitted.
#' @param alternative is a character string, one of \code{"greater"},
#' \code{"less"}, or \code{"two.sided"}, or the initial letter of each,
#' indicating the specification of the alternative hypothesis. For one-sample
#' tests, \code{alternative} refers to the true median of the parent population
#' in relation to the hypothesized value of the mean.
#' @param paired a logical indicating whether you want a paired z-test.
#' @param mu a number specifying the hypothesized mean of the population.
#' @param conf.level confidence level for the interval computation.
#' @param trim the fraction (0 to 0.5) of observations to be trimmed from each
#' end of x before the mean is computed. Values of trim outside that range are
#' taken as the nearest endpoint.
#' @param formula a formula of the form \code{lhs ~ rhs} where \code{lhs} gives
#' the data values and rhs the corresponding groups.
#' @param data an optional matrix or data frame (or similar: see
#' \code{\link{model.frame}}) containing the variables in the formula
#' \code{formula}.  By default the variables are taken from
#' \code{environment(formula)}.
#' @param subset an optional vector specifying a subset of observations to be
#' used.
#' @param na.action a function which indicates what should happen when the data
#' contain NAs. Defaults to \code{getOption("na.action")}.
#' @param \dots further arguments to be passed to or from methods.
#' @return An object of class \code{htest} containing the following components:
#' \item{statistic}{the value of the t-statistic.} \item{parameter}{the degrees
#' of freedom for the t-statistic and the trim percentage used.}
#' \item{p.value}{the p-value for the test.} \item{conf.int}{a confidence
#' interval for the trimmed mean appropriate to the specified alternative
#' hypothesis.} \item{estimate}{the estimated trimmed mean or difference in
#' trimmed means depending on whether it was a one-sample test or a two-sample
#' test. } \item{null.value}{the specified hypothesized value of the trimmed
#' mean or trimmed mean difference depending on whether it was a one-sample
#' test or a two-sample test.} \item{alternative}{a character string describing
#' the alternative hypothesis.} \item{method}{a character string indicating
#' what type of test was performed.} \item{data.name}{a character string giving
#' the name(s) of the data.}
#' 
#' @author Andri Signorell <andri@@signorell.net>, based on R-Core code of
#' \code{\link{t.test}}
#' 
#' @seealso \code{\link{t.test}}, \code{\link{print.htest}}
#' @references Wilcox, R. R. (2005) Introduction to robust estimation and
#' hypothesis testing. \emph{Academic Press}.\cr Yuen, K. K. (1974) The
#' two-sample trimmed t for unequal population variances. \emph{Biometrika},
#' 61, 165-170.
#' 
#' @family topic.parametricTests
#' @concept robust statistics
#' @concept trimmed mean
#' @concept mean comparison
#' 
#' @examples
#' 
#' x <- rnorm(25, 100, 5)
#' yuenTTest(x, mu=99)
#' 
#' # the classic interface
#' with(sleep, yuenTTest(extra[group == 1], extra[group == 2]))
#' 
#' # the formula interface
#' yuenTTest(extra ~ group, data = sleep)
#' 
#' 
#' # Stahel (2002), pp. 186, 196  
#' d.tyres <- data.frame(A=c(44.5,55,52.5,50.2,45.3,46.1,52.1,50.5,50.6,49.2),
#'                       B=c(44.9,54.8,55.6,55.2,55.6,47.7,53,49.1,52.3,50.7))
#' with(d.tyres, yuenTTest(A, B, paired=TRUE))
#' 
#' 
#' d.oxen <- data.frame(ext=c(2.7,2.7,1.1,3.0,1.9,3.0,3.8,3.8,0.3,1.9,1.9),
#'                      int=c(6.5,5.4,8.1,3.5,0.5,3.8,6.8,4.9,9.5,6.2,4.1))
#' with(d.oxen, yuenTTest(int, ext, paired=FALSE))
#' 



#' @rdname yuenTTest
#' @export
yuenTTest <- function (x, ...)
  UseMethod("yuenTTest")


#' @rdname yuenTTest
#' @export
yuenTTest.formula <- function (formula, data, subset, na.action, ...)  {
  
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
  DATA <- setNames(split(mf[[response]], g), c("x", "y"))
  y <- do.call("yuenTTest", c(DATA, list(...)))
  y$data.name <- DNAME
  if (length(y$estimate) == 2L)
    names(y$estimate) <- paste("trimmed mean in group", levels(g))
  y
}


#' @rdname yuenTTest
#' @export
yuenTTest.default <- function (x, y = NULL, alternative = c("two.sided", "less", "greater"),
                               mu = 0, paired = FALSE, conf.level = 0.95, trim = 0.2, ...) {
  
  alternative <- match.arg(alternative)
  if (!missing(mu) && (length(mu) != 1 || is.na(mu)))
    stop("'mu' must be a single number")
  if (!missing(conf.level) && (length(conf.level) != 1 || !is.finite(conf.level) ||
                               conf.level < 0 || conf.level > 1))
    stop("'conf.level' must be a single number between 0 and 1")
  if (!is.null(y)) {
    dname <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
    if (paired)
      xok <- yok <- complete.cases(x, y)
    else {
      yok <- !is.na(y)
      xok <- !is.na(x)
    }
    y <- y[yok]
  }
  else {
    dname <- deparse(substitute(x))
    if (paired)
      stop("'y' is missing for paired test")
    xok <- !is.na(x)
    yok <- NULL
  }
  x <- x[xok]
  
  nx <- length(x)
  mx <- mean(x, trim = trim)
  vx <- var(.winsorize(x, val=quantile(x, probs = c(trim, 1-trim), na.rm=TRUE) ))
  
  if (is.null(y) | paired) {
    if (nx < 2)
      stop("not enough 'x' observations")
    
    df <- nx - 2 * floor(trim * nx) - 1
    
    if(paired){
      my <- mean(y, trim = trim)
      vy <- var(.winsorize(y, val=quantile(y, probs = c(trim, 1-trim), na.rm=TRUE)))
      covxy <- var(.winsorize(x, val=quantile(x, probs = c(trim, 1-trim), na.rm=TRUE)), 
                   .winsorize(y, val=quantile(y, probs = c(trim, 1-trim), na.rm=TRUE)))
      stderr <- sqrt( (nx-1) * (vx + vy - 2 * covxy) / ((df + 1) * df) )
    } else {
      stderr <- sqrt(vx) / ((1 - 2 * trim) * sqrt(nx))
    }
    
    if (stderr < 10 * .Machine$double.eps * abs(mx))
      stop("data are essentially constant")
    
    if(paired){
      method <- "Yuen Paired t-test"
      tstat <- (mx - my - mu) / stderr
      estimate <- setNames(mx - my, "difference of trimmed means")
      
    } else {
      method <- "Yuen One Sample t-test"
      tstat <- (mx - mu)/stderr
      estimate <- setNames(mx, "trimmed mean of x")
    }
    
  }
  else {
    ny <- length(y)
    if (nx < 2)
      stop("not enough 'x' observations")
    if (ny < 2)
      stop("not enough 'y' observations")
    my <- mean(y, trim = trim)
    vy <- var(.winsorize(y, val=quantile(y, probs = c(trim, 1-trim), na.rm=TRUE)))
    method <- "Yuen Two Sample t-test"
    estimate <- c(mx, my)
    names(estimate) <- c("trimmed mean of x", "trimmed mean of y")
    
    dfx <- length(x) - 2 * floor(trim * length(x)) - 1
    dfy <- length(y) - 2 * floor(trim * length(y)) - 1
    
    stderrx <- (length(x) - 1) * vx / ((dfx + 1) * dfx)
    stderry <- (length(y) - 1) * vy / ((dfy + 1) * dfy)
    
    df <- (stderrx + stderry)^2 / (stderrx^2 / dfx + stderry^2 / dfy)
    
    stderr <- sqrt(stderrx + stderry)
    
    if (stderr < 10 * .Machine$double.eps * max(abs(mx), abs(my)))
      stop("data are essentially constant")
    tstat <- (mx - my - mu) / stderr
  }
  if (alternative == "less") {
    pval <- pt(tstat, df)
    cint <- c(-Inf, tstat + qt(conf.level, df))
  }
  else if (alternative == "greater") {
    pval <- pt(tstat, df, lower.tail = FALSE)
    cint <- c(tstat - qt(conf.level, df), Inf)
  }
  else {
    pval <- 2 * pt(-abs(tstat), df)
    alpha <- 1 - conf.level
    cint <- qt(1 - alpha/2, df)
    cint <- tstat + c(-cint, cint)
  }
  cint <- mu + cint * stderr
  names(tstat) <- "t"
  names(df) <- "df"
  names(trim) <- "trim"
  names(mu) <- if (paired || !is.null(y))
    "difference in trimmed means"
  else "trimmed mean"
  attr(cint, "conf.level") <- conf.level
  rval <- list(statistic = tstat, parameter = c(df, trim), p.value = pval,
               conf.int = cint, estimate = estimate, null.value = mu,
               alternative = alternative, method = method, data.name = dname)
  class(rval) <- "htest"
  return(rval)
}




# == internal helper functions ===========================================


.yuenTTestB <- function(x, y, trim = 0, conf.level = 0.95, nboot=599
                        , alternative = c("two.sided", "less", "greater"), mu = 0, na.rm = FALSE){
  
  
  TrimSE <- function(x, trim = 0, na.rm = FALSE) {
    
    #  Estimate the standard error of the gamma trimmed mean
    #  The default amount of trimming is trim = 0.2
    
    if(na.rm) x <- na.omit(x)
    
    winvar <- var(.winsorize(x, val=quantile(x, probs = c(trim, 1-trim), na.rm=TRUE)))
    
    trimse <- sqrt(winvar) / ((1 - 2 * trim) * sqrt(length(x)))
    trimse
  }
  
  
  alternative <- match.arg(alternative)
  method <- "Yuen Two Sample bootstrap t-test"
  dname <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
  
  if(na.rm) x <- na.omit(x)
  if(na.rm) y <- na.omit(y)
  
  meanx <- mean(x, trim = trim)
  meany <- mean(y, trim = trim)
  
  tstat <- (meanx - meany ) / sqrt(TrimSE(x, trim = trim)^2 + TrimSE(y, trim = trim)^2)
  
  sampx <- matrix(sample(x - meanx, size=length(x) * nboot, replace=TRUE), nrow=nboot)
  sampy <- matrix(sample(y - meany, size=length(y) * nboot, replace=TRUE), nrow=nboot)
  
  top <- apply(sampx, 1, mean, trim) - apply(sampy, 1, mean, trim)
  botx <- apply(sampx, 1, TrimSE, trim)
  boty <- apply(sampy, 1, TrimSE, trim)
  tval <- top / sqrt(botx^2 + boty^2)
  
  
  alpha <- 1 - conf.level
  se <- sqrt((TrimSE(x, trim = trim))^2 + (TrimSE(y, trim = trim))^2)
  
  if(alternative == "two.sided") {
    tval <- abs(tval)
    icrit <- floor((1 - alpha) * nboot + .5)
    cint <- meanx - meany + c(-1, 1) * tval[icrit] * se
    pval <- (sum(abs(tstat) <= abs(tval))) / nboot
    
  } else {
    tval <- sort(tval)
    ibot <- floor(alpha/2 * nboot + .5)
    itop <- floor((1 - alpha/2) * nboot + .5)
    cint <- meanx - meany - tval[c(itop, ibot)] * se
    
  }
  
  names(tstat) <- "t"
  names(mu) <- "difference in means"
  estimate <- c(meanx, meany)
  names(estimate) <- c("mean of x", "mean of y")
  
  attr(cint, "conf.level") <- conf.level
  rval <- list(statistic = tstat, p.value = pval,
               conf.int = cint, estimate = estimate, null.value = mu,
               alternative = alternative, method = method, data.name = dname)
  class(rval) <- "htest"
  return(rval)
  
}



# == internal helper functions ===========================================

# this is verbal copy from DescToolsX::winsorize()


.winsorize <- function (x, val = quantile(x, probs = c(0.05, 0.95), 
                                          na.rm = FALSE)) {
  x[x < val[1L]] <- val[1L]
  x[x > val[2L]] <- val[2L]
  return(x)
}


