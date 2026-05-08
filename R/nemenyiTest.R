
#' Nemenyi Test 
#' 
#' A nonparametric post hoc test for multiple pairwise comparisons following 
#' a significant Kruskal-Wallis or Friedman test, based on rank differences.
#' 
#' Performs Nemenyi's test of multiple comparisons. 
#' 
#' 
#' ToDo!! ****************************************************************
#' Tell when to use this test. References needed! Nemenyi proposed a
#' test based on rank sums and the application of the family-wise error method
#' to control Type I error inflation, if multiple comparisons are done. The
#' Tukey and Kramer approach uses mean rank sums and can be employed for
#' equally as well as unequally sized samples without ties.
#' ToDo!! ****************************************************************
#' 
#' @name nemenyiTest
#' @aliases nemenyiTest nemenyiTest.default nemenyiTest.formula
#' @param x a numeric vector of data values, or a list of numeric data vectors.
#' @param g a vector or factor object giving the group for the corresponding
#' elements of \code{x}.  Ignored if \code{x} is a list.
#' @param dist the distribution used for the test. Can be \code{tukey}
#' (default) or \code{chisq}.
#' @param out.list logical, defining if the output should be organized in
#' listform.
#' @param formula a formula of the form \code{lhs ~ rhs} where \code{lhs} gives
#' the data values and \code{rhs} the corresponding groups.
#' @param data an optional matrix or data frame (or similar: see
#' \code{\link{model.frame}}) containing the variables in the formula
#' \code{formula}.  By default the variables are taken from
#' \code{environment(formula)}.
#' @param subset an optional vector specifying a subset of observations to be
#' used.
#' @param na.action a function which indicates what should happen when the data
#' contain \code{NA}s.  Defaults to \code{getOption("na.action")}.
#' @param \dots further arguments to be passed to or from methods.
#' @return A list of class \code{htest}, containing the following components:
#' \item{statistic}{ Nemenyi test} \item{p.value}{ the p-value for the test}
#' \item{null.value}{is the value of the median specified by the null
#' hypothesis. This equals the input argument \code{mu}. } \item{alternative}{a
#' character string describing the alternative hypothesis.} \item{method}{ the
#' type of test applied} \item{data.name}{a character string giving the names
#' of the data.}
#' 
#' @seealso \code{\link{dunnTest}}, \code{\link{conoverTest}} 
#' 
#' @references Nemenyi, P. B. (1963) \emph{Distribution-Free Multiple
#' Comparisons} New York, State University of New York, Downstate Medical
#' Center
#' 
#' Hollander, M., Wolfe, D.A. (1999) \emph{Nonparametric Statistical Methods}
#' New York, Wiley, pp. 787
#' 
#' Friedman, M. (1937) The use of ranks to avoid the assumption of normality
#' implicit in the analysis of variance \emph{Journal of the American
#' Statistical Association}, 32:675-701
#' 
#' Friedman, M. (1940) A comparison of alternative tests of significance for
#' the problem of m rankings \emph{Annals of Mathematical Statistics}, 11:86-92
#' 
#' @examples
#' 
#' ## Hollander & Wolfe (1973), 116.
#' ## Mucociliary efficiency from the rate of removal of dust in normal
#' ##  subjects, subjects with obstructive airway disease, and subjects
#' ##  with asbestosis.
#' x <- c(2.9, 3.0, 2.5, 2.6, 3.2) # normal subjects
#' y <- c(3.8, 2.7, 4.0, 2.4)      # with obstructive airway disease
#' z <- c(2.8, 3.4, 3.7, 2.2, 2.0) # with asbestosis
#' 
#' nemenyiTest(list(x, y, z))
#' 
#' ## Equivalently,
#' x <- c(x, y, z)
#' g <- factor(rep(1:3, c(5, 4, 5)),
#'             labels = c("Normal subjects",
#'                        "Subjects with obstructive airway disease",
#'                        "Subjects with asbestosis"))
#' 
#' nemenyiTest(x, g)
#' 
#' ## Formula interface.
#' boxplot(Ozone ~ Month, data = airquality)
#' nemenyiTest(Ozone ~ Month, data = airquality)
#' 
#' # Hedderich & Sachs, 2012, p. 555
#' d.frm <- data.frame(x=c(28,30,33,35,38,41, 36,39,40,43,45,50, 44,45,47,49,53,54),
#'                     g=c(rep(LETTERS[1:3], each=6)), stringsAsFactors=TRUE)
#' 
#' nemenyiTest(x~g, d.frm)



#' @rdname nemenyiTest
#' @family test.posthoc
#' @concept multiple-comparisons
#' @concept nonparametric
#' @concept hypothesis-testing
#'
#'
#' @export
nemenyiTest <- function (x, ...)
  UseMethod("nemenyiTest")


#' @rdname nemenyiTest
#' @export
nemenyiTest.formula <- function (formula, data, subset, na.action, ...) {
  
  if (missing(formula) || (length(formula) != 3L) || (length(attr(terms(formula[-2L]),
                                                                  "term.labels")) != 1L))
    stop("'formula' missing or incorrect")
  m <- match.call(expand.dots = FALSE)
  if (is.matrix(eval(m$data, parent.frame())))
    m$data <- as.data.frame(data)
  m[[1L]] <- quote(stats::model.frame)
  m$... <- NULL
  mf <- eval(m, parent.frame())
  if (length(mf) > 2L)
    stop("'formula' should be of the form response ~ group")
  DNAME <- paste(names(mf), collapse = " by ")
  names(mf) <- NULL
  response <- attr(attr(mf, "terms"), "response")
  y <- do.call("nemenyiTest", c(as.list(mf), list(...)))
  y$data.name <- DNAME
  y
}




#' @rdname nemenyiTest
#' @export
nemenyiTest.default <- function (x, g,
                                 dist = c("tukey", "chisq"), out.list = TRUE, ...) {
  
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
    g <- factor(rep.int(seq_len(k), l))
    x <- unlist(x)
  }
  else {
    if (length(x) != length(g)) 
      stop("'x' and 'g' must have the same length")
    DNAME <- paste(deparse1(substitute(x)), "and", deparse1(substitute(g)))
    OK <- complete.cases(x, g)
    x <- x[OK]
    g <- g[OK]
    g <- factor(g)
    k <- nlevels(g)
    if (k < 2L) 
      stop("all observations are in the same group")
  }
  N <- length(x)
  if (N < 2L) 
    stop("not enough observations")  
  
  
  dist <- match.arg(dist, c("tukey", "chisq"))
  
  nms <- levels(g)
  
  n <- tapply(g, g, length)
  rnk <- rank(x)
  mrnk <- tapply(rnk, g, mean)
  
  tau <- table(rnk[allDuplicated(rnk)])
  tiesadj <- min(1, 1 - sum(tau^3 - tau) / (N^3 - N))
  mrnkdiff <- outer(mrnk, mrnk, "-")
  
  if(dist == "chisq"){
    chi <- mrnkdiff^2 / ((N*(N+1)/12) * outer(1/n, 1/n, "+"))
    pvals <- pchisq(tiesadj * chi, df=k-1, lower.tail=FALSE)
  } else {
    z <- abs(mrnkdiff) / sqrt( (N*(N+1)/12) * outer(1/n, 1/n, "+"))
    pvals <- ptukey(z * sqrt(2), nmeans=k, df=Inf, lower.tail=FALSE)
  }
  
  
  keep <- lower.tri(pvals)
  pvals <- pvals[keep]
  m <- sum(keep)
  
  out <- list()
  
  # no p.adjustment in this test
  # pvals <- p.adjust(pvals, method=method)
  method.str <- "none" #method
  
  if(out.list){
    dnames <- list(NULL, c("mean rank diff", "pval"))
    if (!is.null(nms))
      dnames[[1L]] <- outer(nms, nms, paste, sep = "-")[keep]
    out[[1]] <- array(c(mrnkdiff[keep], pvals), c(length(mrnkdiff[keep]), 2L), dnames)
    
  } else {
    out[[1]] <- matrix(NA, nrow=length(nms), ncol=length(nms))
    out[[1]][lower.tri(out[[1]], diag = FALSE)] <- pvals
    dimnames(out[[1]]) <- list(nms, nms)
    out[[1]] <- out[[1]][-1, -ncol(out[[1]])]
    
  }
  
  class(out) <- c("DunnTest")
  attr(out, "main") <- gettextf("Nemenyi's test of multiple comparisons for independent samples (%s) ", dist)
  attr(out, "method") <- method.str
  attr(out, "out.list") <- out.list
  
  return(out)
  
}


