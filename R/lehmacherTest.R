
#' Lehmacher's Test for Marginal Homogenity 
#' 
#' A nonparametric test for marginal homogeneity in square contingency tables 
#' for dependent samples, based on a normal approximation of the cell 
#' frequency differences.
#' 
#' Performs Lehmacher's chi-squared test for marginal homogenity in a symmetric
#' two-dimensional contingency table. 
#' 
#' The null is that the probabilities of being classified into cells \verb{[i,j]} and
#' \verb{[j,i]} are the same.
#' 
#' If x is a matrix, it is taken as a two-dimensional contingency table, and
#' hence its entries should be nonnegative integers. Otherwise, both x and y
#' must be vectors or factors of the same length. Incomplete cases are removed,
#' vectors are coerced into factors, and the contingency table is computed from
#' these. 
#' 
#' @name lehmacherTest
#' @aliases lehmacherTest print.mtest
#' @param x either a two-dimensional contingency table in matrix form, or a
#' factor object. 
#' @param y a factor object; ignored if x is a matrix. 
#' @param digits a non-null value for digits specifies the minimum number of
#' significant digits to be printed in values. See details in
#' \code{\link{print.default}}.
#' @param \dots further arguments to be passed to or from other methods. They
#' are ignored in this function.
#' @return A list with class \code{"mtest"} containing the following
#' components: \item{statistic}{a vector with the value of the test
#' statistics.} \item{parameter}{the degrees of freedom, which is always 1 in
#' lehmacherTest.} \item{p.value}{a vector with the p-values of the single
#' tests.} \item{p.value.corr}{a vector with the "hochberg" adjusted p-values
#' of the single tests. (See \code{\link{p.adjust}})} \item{method}{a character
#' string indicating what type of test was performed.} \item{data.name}{a
#' character string giving the name of the data.}
#' 
#' @seealso \code{\link{mcnemar.test}} (resp. BowkerTest for a CxC-matrix),
#' \code{\link{stuartMaxwellTest}}, \code{\link{woolfTest}} 
#' @references Lehmacher, W. (1980) Simultaneous sign tests for marginal
#' homogeneity of square contingency tables \emph{Biometrical Journal}, Volume
#' 22, Issue 8, pages 795-798
#' 
#' @examples
#' 
#' x <- matrix(c(400,40,20,10, 
#'               50,300,60,20, 
#'               10,40,120,5, 
#'               5,90,50,80), nrow=4, byrow=TRUE)
#'               
#' lehmacherTest(x)
#' 
#' 


#' @rdname lehmacherTest
#' @family test.marginal
#' @concept hypothesis-testing
#' @concept table-manipulation
#' @concept nonparametric
#'
#'
#' @export
lehmacherTest <- function(x, y = NULL) {
  
  if (is.matrix(x)) {
    r <- nrow(x)
    if ((r < 2L) || (ncol(x) != r))
      stop("'x' must be square with at least two rows and columns")
    if (any(x < 0, na.rm = TRUE) || any(!is.finite(x)))
      stop("all entries of 'x' must be nonnegative and finite")
    if (any(x != round(x)))
      warning("'x' contains non-integer counts", call. = FALSE)
    DNAME <- deparse(substitute(x))
    
  } else {
    
    if (is.null(y))
      stop("if 'x' is not a matrix, 'y' must be given")
    if (length(x) != length(y))
      stop("'x' and 'y' must have the same length")
    
    DNAME <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
    OK  <- complete.cases(x, y)
    
    ## Unite levels to ensure consistent category set
    lev <- unique(c(as.character(x[OK]), as.character(y[OK])))
    x   <- factor(x[OK], levels = lev)
    y   <- factor(y[OK], levels = lev)
    
    r <- nlevels(x)
    if (r < 2L)
      stop("'x' and 'y' must have at least 2 distinct levels")
    
    x <- table(x, y)
  }
  
  rsum <- rowSums(x)
  csum <- colSums(x)
  
  ## Guard against perfect agreement (denom = 0 -> statistic = 0)
  denom     <- rsum + csum - 2 * diag(x)
  STATISTIC <- ifelse(denom > 0, (rsum - csum)^2 / denom, 0)
  
  PVAL <- pchisq(STATISTIC, df = 1L, lower.tail = FALSE)
  
  structure(
    list(
      statistic    = STATISTIC,
      parameter    = c(df = 1L),
      p.value      = PVAL,
      p.value.corr = p.adjust(PVAL, "hochberg"),
      method       = "Lehmacher-Test on Marginal Homogeneity",
      data.name    = DNAME
    ),
    class = "mtest"
  )
}



#' @rdname lehmacherTest
#' @export
print.mtest <- function(x, digits = 4L, ...) {
  
  cat("\n")
  cat(strwrap(x$method, prefix = "\t"), sep = "\n")
  cat("\n")
  cat("data:  ", x$data.name, "\n", sep = "")
  
  out <- cbind(
    format(round(x$statistic, 4)),
    format.pval(x$p.value,      digits = digits),
    format.pval(x$p.value.corr, digits = digits),
    symnum(x$p.value.corr,
           corr      = FALSE,
           na        = FALSE,
           cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
           symbols   = c("***", "**", "*", ".", " "))
  )
  
  colnames(out) <- c("X-squared", "pval", "pval adj", " ")
  rownames(out) <- if (!is.null(names(x$statistic)))
    names(x$statistic)
  else
    seq_along(x$statistic)
  
  print.default(out, digits = 3L, quote = FALSE, right = TRUE)
  
  cat("\n---\nSignif. codes:  ",
      "0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n\n")
  invisible(x)
}
