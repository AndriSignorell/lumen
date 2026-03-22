
#' Exact Version of Jonckheere-Terpstra Test
#' 
#' A nonparametric test for ordered alternatives across k independent groups, 
#' assessing whether a monotonic trend exists in the location parameter 
#' across ordered groups.
#' 
#' Jonckheere-Terpstra test to test for ordered differences among classes.
#' 
#' jonckheereTerpstraTest is the exact (permutation) version of the
#' Jonckheere-Terpstra test.  It uses the statistic \deqn{\sum_{k<l} \sum_{ij}
#' I(X_{ik} < X_{jl}) + 0.5 I(X_{ik} = X_{jl}),} where \eqn{i, j} are
#' observations in groups \eqn{k} and \eqn{l} respectively.  The asymptotic
#' version is equivalent to \code{cor.test(x, g, method="k")}. The exact
#' calculation requires that there be no ties and that the sample size is less
#' than 100. When data are tied and sample size is at most 100 permutation
#' p-value is returned.\cr
#' 
#' If x is a list, its elements are taken as the samples to be compared, and
#' hence have to be numeric data vectors.  In this case, g is ignored, and one
#' can simply use jonckheereTerpstraTest(x) to perform the test.  If the
#' samples are not yet contained in a list, use jonckheereTerpstraTest(list(x,
#' ...)). \cr
#' 
#' Otherwise, \code{x} must be a numeric data vector, and \code{g} must be a
#' vector or factor object of the same length as \code{x} giving the group for
#' the corresponding elements of \code{x}.
#' 
#' @name jonckheereTerpstraTest
#' @aliases jonckheereTerpstraTest jonckheereTerpstraTest.default
#' jonckheereTerpstraTest.formula
#' @param x a numeric vector of data values, or a list of numeric data vectors.
#' @param g a vector or factor object giving the group for the corresponding
#' elements of x. Ignored if x is a list.
#' @param alternative means are monotonic (\code{two.sided}),
#' \code{increasing}, or \code{decreasing}
#' @param nperm number of permutations for the reference distribution.  The
#' default is \code{NULL} in which case the permutation p-value is not
#' computed. It's recommended to set \code{nperm} to 1000 or higher if
#' permutation p-value is desired.
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
#' @param exact logical, defining if the exact test should be calculated. If
#' left to \code{NULL}, the function uses the exact test up to a samplesize of
#' 100 and falls back to normal approximation for larger samples. The exact
#' procedure can not be applied to samples containing ties.
#' @param \dots further argument to be passed to methods.
#' 
#' @note The function was previously published as \code{jonckheere.test()} in
#' the \pkg{clinfun} package and has been integrated here without logical
#' changes. Some argument checks and a formula interface were added.
#' 
#' @note
#' Original C code by Venkatraman E. Seshan. Rewritten in C++ with an 
#' adapted R interface to conform to package standards.
#'  
#' @references Jonckheere, A. R. (1954). A distribution-free k-sample test
#' again ordered alternatives. \emph{Biometrika} 41:133-145.
#' 
#' Terpstra, T. J. (1952). The asymptotic normality and consistency of
#' Kendall's test against trend, when ties are present in one ranking.
#' \emph{Indagationes Mathematicae} 14:327-333.
#' 
#' @family topic.nonparametricTests
#' @concept ordered alternatives
#' @concept rank-based
#' 
#' @examples
#' 
#' set.seed(1234)
#' g <- ordered(rep(1:5, rep(10,5)))
#' x <- rnorm(50) + 0.3 * as.numeric(g)
#' 
#' jonckheereTerpstraTest(x, g)
#' 
#' x[1:2] <- mean(x[1:2]) # tied data
#' 
#' jonckheereTerpstraTest(x, g)
#' jonckheereTerpstraTest(x, g, nperm=5000)
#' 
#' # Duller, S. 222
#' coffee <- list(
#'   c_4=c(447,396,383,410),
#'   c_2=c(438,521,468,391,504,472),
#'   c_0=c(513,543,506,489,407))  
#' 
#' # the list interface:
#' jonckheereTerpstraTest(coffee)
#' 
#' # the formula interface
#' breaking <- data.frame(
#'   speed=c(20,25,25,25,25,30,30,30,35,35),
#'   distance=c(48,33,59,48,56,60,101,67,85,107))
#' 
#' jonckheereTerpstraTest(distance ~ speed, breaking)
#' 
#' 


#' @rdname jonckheereTerpstraTest
#' @export
jonckheereTerpstraTest <- function (x, ...)  UseMethod("jonckheereTerpstraTest")


#' @rdname jonckheereTerpstraTest
#' @export
jonckheereTerpstraTest.formula <- local({
  
  # super elegant formula implementation
  # in fact we need nothing other, than is already implemented in
  # t.test.formula, besides the last call of zTest() instead of t.test()
  
  tf <- getS3method("kruskal.test", "formula")
  
  new_body <- .replace_text_calls(body(tf), 
                                  old="kruskal.test", 
                                  new="jonckheereTerpstraTest")
  
  new_fun <- tf
  body(new_fun) <- new_body
  
  new_fun
  
})




#' @rdname jonckheereTerpstraTest
#' @export
jonckheereTerpstraTest.default <- function (x, g, 
                                            alternative = c("two.sided", "increasing", "decreasing"), 
                                            nperm=NULL, exact=NULL,...) {
  
  
  if (is.list(x)) {
    if (length(x) < 2L) 
      stop("'x' must be a list with at least 2 elements")
    if (!missing(g)) 
      warning("'x' is a list, so ignoring argument 'g'")
    DNAME <- deparse1(substitute(x))
    x <- lapply(x, function(u) u <- u[complete.cases(u)])
    if (!all(sapply(x, is.numeric))) 
      warning("some elements of 'x' are not numeric and will be coerced to numeric")
    k <- length(x)
    l <- lengths(x)
    if (any(l == 0L)) 
      stop("all groups must contain data")
    g <- ordered(rep.int(seq_len(k), l))
    x <- unlist(x)
    
  } else {
    
    if (length(x) != length(g)) 
      stop("'x' and 'g' must have the same length")
    DNAME <- paste(deparse1(substitute(x)), "by", 
                   deparse1(substitute(g)))
    OK <- complete.cases(x, g)
    x <- x[OK]
    g <- g[OK]
    g <- ordered(g)
    k <- nlevels(g)
    if (k < 2L) 
      stop("all observations are in the same group")
  }
  n <- length(x)
  if (n < 2L) 
    stop("not enough observations")
  
  
  # start calculating
  
  jtpdf <- function(gsize) {
    gsize  <- as.integer(gsize)
    ng     <- length(gsize)
    cgsize <- rev(cumsum(rev(gsize)))
    mxsum  <- sum(gsize[-ng] * cgsize[-1]) + 1
    
    jtpdf_cpp(
      mxsum,
      cgsize,
      numeric(mxsum),
      numeric(mxsum)
    )
  }
  
  
  if(!is.numeric(g) & !is.ordered(g)) stop("group should be numeric or ordered factor")
  
  alternative <- match.arg(alternative)
  
  
  jtperm.p <- function(x, ng, gsize, cgsize, alternative, nperm) {
    # this function computes the pdf using the convolution by Mark van de Wiel
    
    n <- length(x)
    pjtrsum <- rep(0, nperm)
    for (np in 1:nperm){
      jtrsum <- 0
      for(i in 1L:(ng-1)) {
        na <- gsize[i]
        nb <- n-cgsize[i+1]
        # this jtrsum will be small if data are increasing and large if decreasing
        jtrsum <- jtrsum + sum(rank(x[(cgsize[i]+1):n])[1:na]) - na*(na+1)/2
      }
      pjtrsum[np] <- jtrsum
      # permute the data; this way the first value is the original statistic
      x <- sample(x)
    }
    # one-sided p-values
    # number of permuted values at least as small as original
    iPVAL <- sum(pjtrsum <= pjtrsum[1])/nperm
    # number of permuted values at least as large as original
    dPVAL <- sum(pjtrsum >= pjtrsum[1])/nperm
    # return p-value for the alternative of interest
    PVAL <- switch(alternative,
                   "two.sided" = 2*min(iPVAL, dPVAL, 1),
                   "increasing" = iPVAL,
                   "decreasing" = dPVAL)
    PVAL
  }
  
  
  # Alternative for the JT-Statistic
  # JT <- function(z){
  #
  #   w <- function(x, y){
  #     # verbatim from wilcox.test STATISTIC
  #     r <- rank(c(x, y))
  #     n.x <- as.double(length(x))
  #     n.y <- as.double(length(y))
  #     STATISTIC <- c(W = sum(r[seq_along(x)]) - n.x * (n.x + 1)/2)
  #     STATISTIC
  #   }
  #
  #   k <- length(z)
  #   u <- 0
  #
  #   for(i in 2:k){
  #     for(j in 1:(i-1))	{
  #       u <- u + w(z[[i]], z[[j]])
  #     } }
  #   u
  # }
  
  
  # see:
  #   library(NSM3)
  #   HOllander Wolfe pp 218
  #   piece <- c(40,35,38,43,44,41, 38,40,47,44,40,42, 48,40,45,43,46,44)
  #   grp <- factor(rep(c("ctrl","A","B"), each=6), ordered=TRUE, levels=c("ctrl","A","B"))
  #
  #   jonckheereTerpstraTest(piece, grp)
  #   pJCK(piece, grp)
  
  
  METHOD <- "Jonckheere-Terpstra test"
  PERM <- !missing(nperm)
  n <- length(x)
  if(length(g) != n) stop("lengths of data values and group don't match")
  TIES <- length(unique(x)) != n
  gsize <- table(g)
  ng <- length(gsize)
  cgsize <- c(0,cumsum(gsize))
  x <- x[order(g)]
  jtrsum <- jtmean <- jtvar <- 0
  for(i in 1:(ng-1)) {
    na <- gsize[i]
    nb <- n-cgsize[i+1]
    jtrsum <- jtrsum + sum(rank(x[(cgsize[i]+1):n])[1:na]) - na*(na+1)/2
    jtmean <- jtmean + na*nb/2
    jtvar <- jtvar + na*nb*(na+nb+1)/12
  }
  # this jtrsum will be small if data are increasing and large if decreasing
  # to reverse this use 2*jtmean - jtrsum
  jtrsum <- 2*jtmean - jtrsum
  STATISTIC <- jtrsum
  names(STATISTIC) <- "JT"
  
  if(is.null(exact)) {
    exact <- !(n > 100 | TIES)
    if(!exact)
      warning("Sample size > 100 or data with ties \n p-value based on normal approximation. Specify nperm for permutation p-value")
  }
  
  if(exact & TIES){
    warning("Sample data with ties \n p-value based on normal approximation. Specify nperm for permutation p-value.")
    exact <- FALSE
  }
  
  if (PERM) {
    PVAL <- jtperm.p(x, ng, gsize, cgsize, alternative, nperm)
    
  } else {
    if(!exact){
      zstat <- (STATISTIC - jtmean) / sqrt(jtvar)
      PVAL <- pnorm(zstat)
      PVAL <- switch(alternative,
                     "two.sided" = 2 * min(PVAL, 1-PVAL, 1),
                     "increasing" = 1-PVAL,
                     "decreasing" = PVAL)
      
    } else {
      dPVAL <- sum(jtpdf(gsize)[1:(jtrsum+1)])
      iPVAL <- 1-sum(jtpdf(gsize)[1:(jtrsum)])
      PVAL <- switch(alternative,
                     "two.sided" = 2 * min(iPVAL, dPVAL, 1),
                     "increasing" = iPVAL,
                     "decreasing" = dPVAL)
      
    }
  }
  
  RVAL <- list(statistic = STATISTIC,
               p.value = as.numeric(PVAL),
               alternative = alternative,
               method = METHOD,
               data.name = DNAME)
  class(RVAL) <- "htest"
  RVAL
  
}



# jonckheereTerpstraTest.formula <- function (formula, data, subset, na.action, ...) {
#   
#   if (missing(formula) || (length(formula) != 3L)) 
#     stop("'formula' missing or incorrect")
#   m <- match.call(expand.dots = FALSE)
#   if (is.matrix(eval(m$data, parent.frame()))) 
#     m$data <- as.data.frame(data)
#   m[[1L]] <- quote(stats::model.frame)
#   m$... <- NULL
#   mf <- eval(m, parent.frame())
#   if (length(mf) > 2L) 
#     stop("'formula' should be of the form response ~ group")
#   DNAME <- paste(names(mf), collapse = " by ")
#   names(mf) <- NULL
#   y <- do.call("jonckheereTerpstraTest", c(as.list(mf), list(...)))
#   y$data.name <- DNAME
#   y
#   
# }
# 
# 
