
#' Bhapkar Marginal Homogeneity Test 
#' 
#' A nonparametric test for marginal homogeneity in square contingency 
#' tables for dependent samples, similar to the Stuart-Maxwell test but 
#' based on a different test statistic and generally slightly more powerful.
#' 
#' Bhapkar’s test (Bhapkar, 1966) is used to assess marginal homogeneity 
#' in square contingency tables. It is based on the asymptotic 
#' normality of marginal proportions and is closely related to 
#' the generalized McNemar test, as implemented in \code{stuartMaxwellTest}. 
#' The two tests differ primarily in the construction of the 
#' variance–covariance matrix.
#' 
#' The Bhapkar and Stuart–Maxwell tests are asymptotically 
#' equivalent (Keefe, 1982), meaning that for large sample 
#' sizes they yield the same chi-square statistic. For finite 
#' samples, however, the Bhapkar test is generally more powerful 
#' and is therefore preferred in practice. In particular, the difference 
#' between the two tests lies in the estimation of the variance–covariance 
#' matrix of the marginal proportions.
#' 
#' @name bhapkarTest
#' @docType data
#' @param x either a 2-way \eqn{k \times k}{k x k} contingency table in matrix
#' form, or a factor. 
#' @param y a factor with the same levels as \code{x}; ignored if \code{x} is a
#' matrix. 
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
#' @examples
#' 
#' # Source: https://john-uebersax.com-us.com/stat/mcnemar.htm#bhapkar
#' mc <- as.table(matrix(c(20,3,0,10,30,5,5,15,40), nrow=3))
#' 
#' bhapkarTest(mc)
#' 



#' @rdname bhapkarTest
#' @family test.marginal
#' @concept hypothesis-testing
#' @concept table-manipulation
#' @concept nonparametric
#'
#'
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

