
#' Bhapkar Marginal Homogeneity Test 
#' 
#' Bhapkar (1966) tested marginal homogeneity by exploiting the asymptotic
#' normality of marginal proportion, and so this test is also called Bhapkar's
#' test. The idea of constructing test statistic is similar to the one of
#' generalized McNemar's test statistic used in
#' \code{\link{stuartMaxwellTest}}, and the major difference lies in the
#' calculation of elements in variance-covariance matrix.  
#' 
#' Although the Bhapkar and Stuart-Maxwell tests are asymptotically equivalent
#' (Keefe, 1982). Generally, the Bhapkar (1966) test is a more powerful
#' alternative to the Stuart-Maxwell test. With a large N, both will produce
#' the same Chi-square value. As the Bhapkar test is more powerful, it is
#' preferred.
#' 
#' @name bhapkarTest
#' @docType data
#' @param x either a 2-way \eqn{k \times k}{k x k} contingency table in matrix
#' form, or a factor. 
#' @param y a factor with the same levels as \code{x}; ignored if \code{x} is a
#' matrix. 
#' @author Andri Signorell <andri@@signorell.net> 
#' 
#' @seealso \code{\link{mcnemar.test}},
#' \code{\link{chisq.test}}, \code{\link{mhChisqTest}},
#' \code{\link{breslowDayTest}} 
#' 
#' @references Bhapkar V.P. (1966) A note on the equivalence of two test
#' criteria for hypotheses in categorical data. \emph{Journal of the American
#' Statistical Association}, 61: 228-235.
#' 
#' Ireland C.T., Ku H.H., and Kullback S. (1969) Symmetry and marginal
#' homogeneity of an r x r contingency table. \emph{Journal of the American
#' Statistical Association}, 64: 1323-1341. 
#' 
#' Keefe T.J. (1982) On the relationship between two tests for homogeneity of
#' the marginal distributions in a two-way classification. \emph{Biometrika},
#' 69: 683-684.
#' 
#' Sun X., Yang Z. (2008) Generalized McNemar's Test for Homogeneity of the
#' Marginal Distributions. \emph{SAS Global Forum 2008: Statistics and Data
#' Analysis}, Paper 382-208.
#' 
#' @family topic.contingencyTests
#' @concept marginal homogeneity
#' 
#' @examples
#' 
#' # Source: https://john-uebersax.com-us.com/stat/mcnemar.htm#bhapkar
#' mc <- as.table(matrix(c(20,3,0,10,30,5,5,15,40), nrow=3))
#' 
#' bhapkarTest(mc)
#' 



#' @rdname bhapkarTest
#' @export
bhapkarTest <- function(x, y = NULL){
  
  # https://support.sas.com/resources/papers/proceedings/pdfs/sgf2008/382-2008.pdf
  
  if (is.matrix(x)) {
    DNAME <- deparse(substitute(x))
  } else {
    DNAME <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
  }
  
  
  res <- stuartMaxwellTest(x=x, y=y)
  
  STATISTIC <- res$statistic
  PARAMETER <- res$parameter
  
  # new statistic by Bhapkar
  STATISTIC <- STATISTIC/(1-STATISTIC/res$n)
  
  res$statistic <- STATISTIC
  res$p.value <- pchisq(STATISTIC, PARAMETER, lower.tail = FALSE)
  
  res$data.name <- DNAME
  res$method <- "Bhapkar test"
  
  return(res)
  
}

