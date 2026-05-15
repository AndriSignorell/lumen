
#' Stuart-Maxwell Marginal Homogeneity Test 
#' 
#' A nonparametric test for marginal homogeneity in square contingency tables 
#' for dependent samples, generalizing the McNemar test to more than two 
#' categories.
#' 
#' This function computes the marginal homogeneity test for a \eqn{k \times
#' k}{k x k} matrix of assignments of objects to \code{k} categories or two
#' vectors \code{x}, \code{y} of category scores for \code{n} data objects by
#' two raters. The statistic is distributed as \eqn{\chi^2}{chi-square} with
#' \code{k-1} degrees of freedom. \cr It can be viewed as an extension of the
#' McNemar test to \eqn{k \times k}{k x k} table.  
#' 
#' The null is that the probabilities of being classified into cells \verb{[i, j]} 
#' and \verb{[j, i]} are the same.
#' 
#' If \code{x} is a matrix, it is taken as a two-dimensional contingency table,
#' and hence its entries should be nonnegative integers. Otherwise, both x and
#' y must be vectors or factors of the same length and with the same levels.
#' \cr Incomplete cases are removed, vectors are coerced into factors, and the
#' contingency table is computed from these.
#' 
#' If there is perfect agreement for any category k, that category must be
#' omitted in order to invert matrix S.
#' 
#' If for any category \code{k}, all frequencies in row \code{k} and column
#' \code{k} are 0, except possibly for the main diagonal element (e.g., for
#' perfect agreement for category \code{k}, in such cases also the
#' corresponding row and column marginal frequencies would be equal), then the
#' category is not included in the test and should be ignored, say the
#' Stuart-Maxwell test is performed with respect to the remaining categories
#' only. The degree of freedom \code{df} in this case can still be considered
#' \code{k - 1}, where \code{k} is the number of original categories; this
#' treats omitted categories as if they were included but contributed 0 to the
#' value of \eqn{\chi^2}{Chi-square} - a reasonable view since such categories
#' have equal row and column marginals. (See:
#' \url{https://www.john-uebersax.com/stat/mcnemar.htm#stuart})
#' 
#' @name stuartMaxwellTest
#' @param x either a 2-way \eqn{k \times k}{k x k} contingency table in matrix
#' form, or a factor. 
#' @param y a factor with the same levels as x; ignored if x is a matrix. 
#' 
#' @return A list with class \code{"htest"} containing the following
#' components: \item{statistic}{the value of the test statistic.}
#' \item{parameter}{the degrees of freedom.} \item{p.value}{the p-value of the
#' test.} \item{method}{a character string indicating what type of test was
#' performed.} \item{data.name}{a character string giving the name of the
#' data.}
#' @note based on code from Jim Lemon
#' 
#' @seealso \code{\link{bhapkarTest}} for a more powerful alternative to the
#' Stuart-Maxwell test
#' 
#' \code{\link{mcnemar.test}}, \code{\link{chisq.test}},
#' \code{\link{mantelTrendTest}}, \code{\link{breslowDayTest}}
#' 
#' @references Stuart, A (1955) A test for homogeneity of the marginal
#' distributions in a two-way classification. \emph{Biometrika}, 42, 412-416.
#' 
#' Maxwell, A.E. (1970) Comparing the classification of subjects by two
#' independent judges. \emph{British Journal of Psychiatry}, 116, 651-655.
#' 
#' Agresti, A. (2002) \emph{Categorical Data Analysis}. John Wiley & Sons, pp
#' 86 ff.
#' 
#' @examples
#' 
#' # Source: https://john-uebersax.com-us.com/stat/mcnemar.htm#stuart
#' hyp <- as.table(matrix(c(20,3,0,10,30,5,5,15,40), nrow=3))
#' stuartMaxwellTest(hyp)
#' 
#' # same as defined with two vectors
#' d.hyp <- expand.grid(c("A","B","C"), c("A","B","C"))[
#'                      rep(1:9, times = hyp), ]
#' row.names(d.hyp) <- NULL
#' 
#' stuartMaxwellTest(x=d.hyp[,1], y=d.hyp[,2])
#' 
#' 
#' mc <- as.table(matrix(c(
#'          732, 1524, 1575, 1577, 1602, 837, 1554, 1437, 
#'          1672, 1600, 841, 1363, 1385, 1484, 1524, 791), nrow=4))
#' 
#' stuartMaxwellTest(mc)
#' 


#' @rdname stuartMaxwellTest
#' @family test.marginal
#' @concept hypothesis-testing
#' @concept table-manipulation
#' @concept nonparametric
#'
#'


#' @export
stuartMaxwellTest <- function(x, y = NULL) {
  
  if (is.matrix(x)) {
    
    if (!is.numeric(x))
      stop("'x' must be a numeric matrix")
    r <- nrow(x)
    if ((r < 2L) || (ncol(x) != r))
      stop("'x' must be square with at least two rows and columns")
    if (any(x < 0, na.rm = TRUE) || any(!is.finite(x)))
      stop("all entries of 'x' must be nonnegative and finite")
    if (any(x != round(x)))
      warning("'x' contains non-integer counts", call. = FALSE)
    if (!is.null(dimnames(x)) &&
        !identical(dimnames(x)[[1L]], dimnames(x)[[2L]]))
      warning(
        "row and column names differ; ensure rows and columns ",
        "represent the same categories",
        call. = FALSE
      )
    
    DNAME <- deparse(substitute(x))
    
  } else {
    
    if (is.null(y))
      stop("if 'x' is not a matrix, 'y' must be given")
    if (length(x) != length(y))
      stop("'x' and 'y' must have the same length")
    
    DNAME <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
    
    OK  <- complete.cases(x, y)
    
    ## Preserve original data order for levels; treat categories as nominal
    lev <- unique(c(as.character(x[OK]), as.character(y[OK])))
    x   <- factor(x[OK], levels = lev)
    y   <- factor(y[OK], levels = lev)
    
    r <- nlevels(x)
    if (r < 2L)
      stop("'x' and 'y' must have at least 2 distinct levels")
    
    x <- table(x, y)
  }
  
  ## Save original k for df (see Details)
  k_original <- nrow(x)
  
  rowsums <- rowSums(x)
  colsums <- colSums(x)
  
  ## Remove categories with perfect agreement
  equalsums <- diag(x) == rowsums & diag(x) == colsums
  if (any(equalsums)) {
    x <- x[!equalsums, !equalsums, drop = FALSE]
    if (nrow(x) < 2L)
      stop("too many categories with perfect agreement; cannot compute")
    rowsums <- rowSums(x)
    colsums <- colSums(x)
  }
  
  k1  <- nrow(x) - 1L
  smd <- (rowsums - colsums)[seq_len(k1)]
  
  ## Variance-covariance matrix S (Stuart 1955, eq. 3)
  smS <- -(x[seq_len(k1), seq_len(k1)] +
             t(x[seq_len(k1), seq_len(k1)]))
  diag(smS) <- rowsums[seq_len(k1)] +
    colsums[seq_len(k1)] -
    2 * diag(x)[seq_len(k1)]
  
  STATISTIC <- tryCatch(
    drop(smd %*% qr.solve(smS, smd, tol = 1e-7)),
    error = function(e)
      stop(
        "variance-covariance matrix is singular; ",
        "cannot compute statistic",
        call. = FALSE
      )
  )
  
  ## df = k_original - 1: omitted categories (perfect agreement) are
  ## treated as contributing 0 to chi-squared — Stuart-Maxwell convention
  PARAMETER <- c(df = k_original - 1L)
  
  structure(
    list(
      statistic = c("chi-squared" = STATISTIC),
      parameter = PARAMETER,
      p.value   = pchisq(STATISTIC, df = PARAMETER, lower.tail = FALSE),
      method    = "Stuart-Maxwell test for marginal homogeneity",
      data.name = DNAME,
      n         = sum(rowsums)
    ),
    class = "htest"
  )
}


