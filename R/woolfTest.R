
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
#' @family test.contingency
#' @concept hypothesis-testing
#' @concept table-manipulation
#'
#'


#' @export
woolfTest <- function(x) {
  
  DNAME <- deparse(substitute(x))
  
  if (!is.array(x) || length(dim(x)) != 3L || any(dim(x)[1:2] != 2L))
    stop("'x' must be a 2x2xK array")
  
  if (any(x < 0, na.rm = TRUE) || any(!is.finite(x)))
    stop("all entries of 'x' must be nonnegative and finite")
  
  if (any(x != round(x)))
    warning("'x' contains non-integer counts", call. = FALSE)
  
  if (any(x == 0)) {
    message("zero cell count detected; adding 0.5 to all cells")
    x <- x + 0.5
  }
  
  k  <- dim(x)[3L]
  or <- apply(x, 3L, function(x) (x[1,1] * x[2,2]) / (x[1,2] * x[2,1]))
  w  <- apply(x, 3L, function(x) 1 / sum(1 / x))
  o  <- log(or)
  e  <- weighted.mean(o, w)
  
  STATISTIC <- sum(w * (o - e)^2)
  PARAMETER <- k - 1L
  
  structure(
    list(
      statistic = c("X-squared" = STATISTIC),
      parameter = c(df = PARAMETER),
      p.value   = pchisq(STATISTIC, PARAMETER, lower.tail = FALSE),
      method    = "Woolf Test on Homogeneity of Odds Ratios (no 3-Way assoc.)",
      data.name = DNAME,
      observed  = o,
      expected  = e
    ),
    class = "htest"
  )
}


