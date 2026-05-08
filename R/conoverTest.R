
#' Conover's Test of Multiple Comparisons
#' 
#' A nonparametric post hoc test for multiple pairwise comparisons following 
#' a significant Kruskal-Wallis test, based on rank data.
#' 
#' Perform Conover's test of multiple comparisons using rank sums as post hoc
#' test following a significant \code{\link{kruskal.test}}.
#' 
#' \code{conoverTest} performs the post hoc pairwise multiple comparisons
#' procedure appropriate to follow the rejection of a Kruskal-Wallis test.
#' Conover's test is more powerful than Dunn's post hoc multiple comparisons
#' test (\code{\link{dunnTest}}). The interpretation of stochastic dominance
#' requires an assumption that the CDF of one group does not cross the CDF of
#' the other. \cr conoverTest makes m = k(k-1)/2 multiple pairwise comparisons
#' based on the Conover-Iman t-test-statistic for the rank-sum differences:
#' \deqn{\left | \bar{R}_{i}-\bar{R}_{j} \right | > t_{1-\alpha/2, n-k} \cdot
#' \sqrt{ s^2 \cdot \left [ \frac{n-1-\hat{H}^*}{n-k} \right ] \cdot \left [
#' \frac{1}{n_i} + \frac{1}{n_j} \right ] } } with the (tie corrected)
#' statistic of the Kruskal Wallis test \deqn{\hat{H}^* = \frac{\frac{12}{n
#' \cdot (n+1)} \cdot \sum_{i=1}^{k}\frac{R_{i}^2}{n_i} - 3\cdot(n+1) }
#' {1-\frac{\sum_{i=1}^{r} \left ( t_i^3-t_i \right )}{n^3-n}} } and the
#' \eqn{s^2} being \deqn{s^2 = \frac{1}{n-1} \cdot \left [ \sum{R_i^2} - n
#' \cdot \frac{(n+1)^2}{4} \right ]}
#' 
#' If \code{x} is a list, its elements are taken as the samples to be compared,
#' and hence have to be numeric data vectors.  In this case, \code{g} is
#' ignored, and one can simply use \code{conoverTest(x)} to perform the test.
#' If the samples are not yet contained in a list, use
#' \code{conoverTest(list(x, ...))}.
#' 
#' Otherwise, \code{x} must be a numeric data vector, and \code{g} must be a
#' vector or factor object of the same length as \code{x} giving the group for
#' the corresponding elements of \code{x}.
#' 
#' @name conoverTest
#' 
#' @aliases conoverTest conoverTest.default conoverTest.formula
#' @param x a numeric vector of data values, or a list of numeric data vectors.
#' @param g a vector or factor object giving the group for the corresponding
#' elements of \code{x}.  Ignored if \code{x} is a list.
#' @param method the method for adjusting p-values for multiple comparisons.
#' The function is calling \code{\link{p.adjust}} and this parameter is
#' directly passed through.
#' @param alternative a character string specifying the alternative hypothesis,
#' must be one of \code{"two.sided"} (default), \code{"greater"} or
#' \code{"less"}. You can specify just the initial letter.
#' @param out.list logical, indicating if the results should be printed in list
#' mode or as a square matrix. Default is list (TRUE).
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
#' @return A list with class \code{"DunnTest"} containing the following
#' components: \item{res}{an array containing the mean rank differencens and
#' the according p-values}
#' 
#' @seealso \code{\link{dunnTest}}, \code{\link{nemenyiTest}},
#' \code{\link{kruskal.test}}, \code{\link{wilcox.test}},
#' \code{\link{p.adjust}}
#' 
#' @references Conover W. J., Iman R. L. (1979) On multiple-comparisons
#' procedures, \emph{Tech. Rep.} LA-7677-MS, Los Alamos Scientific Laboratory.
#' 
#' Conover, W. J. (1999) Practical Nonparametric Statistics \emph{Wiley},
#' Hoboken, NJ. 3rd edition.
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
#' conoverTest(list(x, y, z))
#' 
#' ## Equivalently,
#' x <- c(x, y, z)
#' g <- factor(rep(1:3, c(5, 4, 5)),
#'             labels = c("Normal subjects",
#'                        "Subjects with obstructive airway disease",
#'                        "Subjects with asbestosis"))
#' 
#' # do the kruskal.test first
#' kruskal.test(x, g)
#' 
#' # ...and the pairwise test afterwards
#' conoverTest(x, g)
#' 
#' ## Formula interface.
#' boxplot(Ozone ~ Month, data = airquality)
#' conoverTest(Ozone ~ Month, data = airquality)
#' 

#' @rdname conoverTest
#' @family test.posthoc
#' @concept multiple-comparisons
#' @concept nonparametric
#' @concept hypothesis-testing
#'
#'
#' @export
conoverTest <- function (x, ...)
  UseMethod("conoverTest")


#' @rdname conoverTest
#' @export
conoverTest.default <- function (x, g,
                                 method = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"),
                                 alternative = c("two.sided", "less", "greater"), out.list = TRUE, ...) {
  
  alternative <- match.arg(alternative)
  
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
  
  
  method <- match.arg(method)
  nms <- levels(g)
  n <- tapply(g, g, length)
  rnk <- rank(x)
  mrnk <- tapply(rnk, g, mean)
  tau <- table(rnk[allDuplicated(rnk)])
  tiesadj <- 1-sum(tau^3 - tau)/(N^3 - N)
  mrnkdiff <- outer(mrnk, mrnk, "-")
  
  # Kruskal-Wallis H statistic
  H <- (12 / (N * (N + 1))) * sum(tapply(rnk, g, sum)^2 / n) - 3 * (N + 1)
  if (tiesadj == 1) {
    s2 <- N * (N + 1) / 12
  } else {
    s2 <-   ( 1 / (N - 1)) * (sum(rnk^2) - (N * (((N + 1)^2) / 4)))
  }
  
  tval <- mrnkdiff/sqrt(s2 * ((N - 1 - H/tiesadj) / (N - k)) * outer(1/n, 1/n, "+"))
  
  if (alternative == "less") {
    pvals <- pt(abs(tval), df=N - k)
    
  } else if (alternative == "greater") {
    pvals <- pt(abs(tval), df=N - k, lower.tail = FALSE)
    
  } else {
    pvals <- 2 * pt(abs(tval), df=N - k, lower.tail = FALSE)
    
  }
  
  keep <- lower.tri(pvals)
  pvals <- pvals[keep]
  m <- sum(keep)
  out <- list()
  pvals <- p.adjust(pvals, method = method)
  method.str <- method
  if (out.list) {
    dnames <- list(NULL, c("mean rank diff", "pval"))
    if (!is.null(nms))
      dnames[[1L]] <- outer(nms, nms, paste, sep = "-")[keep]
    out[[1]] <- array(c(mrnkdiff[keep], pvals), c(length(mrnkdiff[keep]),
                                                  2L), dnames)
  } else {
    out[[1]] <- matrix(NA, nrow = length(nms), ncol = length(nms))
    out[[1]][lower.tri(out[[1]], diag = FALSE)] <- pvals
    dimnames(out[[1]]) <- list(nms, nms)
    out[[1]] <- out[[1]][-1, -ncol(out[[1]])]
  }
  
  class(out) <- c("conoverTest", "DunnTest")
  attr(out, "main") <- gettextf("Conover's test of multiple comparisons : %s ",
                                method.str)
  attr(out, "method") <- method.str
  attr(out, "out.list") <- out.list
  
  return(out)
  
}




#' @rdname conoverTest
#' @export
conoverTest.formula <- function (formula, data, subset, na.action, ...) {
  
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
  y <- do.call("conoverTest", c(as.list(mf), list(...)))
  y$data.name <- DNAME
  
  y
  
}

