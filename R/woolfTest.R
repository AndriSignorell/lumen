
#' Woolf Test For Homogeneity in 2x2xk Tables
#' 
#' A test for homogeneity of odds ratios across several 2×2 contingency 
#' tables, similar to the Breslow-Day test but based on a different 
#' test statistic.
#' 
#' Test for homogeneity on \eqn{2 \times 2 \times k}{2 x 2 x k} tables over
#' strata (i.e., whether the log odds ratios are the same in all strata).
#' 
#' 
#' @name woolfTest
#' @param x a \eqn{2 \times 2 \times k}{2 x 2 x k} table, where the last
#' dimension refers to the strata.
#' @return A list of class \code{"htest"} containing the following components:
#' \item{statistic}{the chi-squared test statistic.} \item{parameter}{degrees
#' of freedom of the approximate chi-squared distribution of the test
#' statistic.} \item{p.value}{\eqn{p}-value for the test.} \item{method}{a
#' character string indicating the type of test performed.} \item{data.name}{a
#' character string giving the name(s) of the data.} \item{observed}{the
#' observed counts.} \item{expected}{the expected counts under the null
#' hypothesis.}
#' @note This function was previously published as \code{woolf_test()} in the
#' \pkg{vcd} package and has been integrated here without logical changes.
#' @author David Meyer, Achim Zeileis, Kurt Hornik, Michael Friendly
#' 
#' @seealso \code{\link{mantelhaen.test}}, \code{\link{breslowDayTest}}
#' 
#' @references Woolf, B. 1955: On estimating the relation between blood group
#' and disease. \emph{Ann. Human Genet.} (London) \bold{19}, 251-253.
#' 
#' @family topic.contingencyTests
#' @concept odds ratio
#' 
#' @examples
#' 
#' migraine <- xtabs(freq ~ .,
#'             cbind(expand.grid(treatment=c("active","placebo"),
#'                                response=c("better","same"),
#'                                gender=c("female","male")),
#'                   freq=c(16,5,11,20,12,7,16,19))
#'             )
#' 
#' woolfTest(migraine)
#' 
#' 
# the VCD package (available via CRAN) has a function called woolf_test()


#' @rdname woolfTest
#' @export
woolfTest <- function(x) {
  
  DNAME <- deparse(substitute(x))
  if (any(x == 0))
    x <- x + 1 / 2
  
  k <- dim(x)[3]
  or <- apply(x, 3, function(x) (x[1,1] * x[2,2]) / (x[1,2] * x[2,1]))
  w <-  apply(x, 3, function(x) 1 / sum(1 / x))
  o <- log(or)
  e <- weighted.mean(log(or), w)
  
  STATISTIC <- sum(w * (o - e)^2)
  PARAMETER <- k - 1
  PVAL <- 1 - pchisq(STATISTIC, PARAMETER)
  METHOD <- "Woolf Test on Homogeneity of Odds Ratios (no 3-Way assoc.)"
  names(STATISTIC) <- "X-squared"
  names(PARAMETER) <- "df"
  structure(list(statistic = STATISTIC, parameter = PARAMETER,
                 p.value = PVAL, method = METHOD, data.name = DNAME, 
                 observed = o, expected = e), 
            class = "htest")
  
}


#' 