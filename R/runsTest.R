

#' Runs Test for Randomness
#' 
#' A nonparametric test for randomness of a sequence, based on the number 
#' of runs (consecutive sequences of identical values or values 
#' above/below a threshold).
#' 
#' Performs a test whether the elements of \code{x} are serially independent -
#' say, whether they occur in a random order - by counting how many runs there
#' are above and below a threshold. If \code{y} is supplied a two sample
#' Wald-Wolfowitz-Test testing the equality of two distributions against
#' general alternatives will be computed.
#' 
#' \bold{The runs test for randomness} \verb{ } is used to test the hypothesis
#' that a series of numbers is random. \cr
#' 
#' For a categorical variable, the number of runs correspond to the number of
#' times the category changes, that is, where \eqn{x_{i}}{x_i} belongs to one
#' category and \eqn{x_{i+1}}{x_(i+1)} belongs to the other. The number of runs
#' is the number of sign changes plus one.\cr
#' 
#' For a numeric variable x containing more than two values, a run is a set of
#' sequential values that are either all above or below a specified cutpoint,
#' typically the median. This is not necessarily the best choice. If another
#' threshold should be used use a code like: \code{runsTest(x > mean(x))}.
#' 
#' The exact distribution of runs and the p-value based on it are described in
#' the manual of SPSS "Exact tests"
#' \url{https://www.sussex.ac.uk/its/pdfs/SPSS_Exact_Tests_21.pdf}.
#' 
#' The normal approximation of the runs test is calculated with the expected
#' number of runs under the null \deqn{\mu_r=\frac{2 n_0 n_1}{n_0 + n_1} + 1}
#' and its variance \deqn{\sigma_r^2 = \frac{2 n_0 n_1 (2 n_0 n_1 - n_0 - n_1)
#' }{(n_0 + n_1)^2 \cdot (n_0 + n_1 - 1)}} as \deqn{\hat{z}=\frac{r - \mu_r +
#' c}{\sigma_r}} where \eqn{n_0, n_1} the number of values below/above the
#' threshold and \eqn{r} the number of runs.
#' 
#' Setting the continuity correction \code{correct = TRUE} will yield the
#' normal approximation as SAS (and SPSS if n < 50) does it, see
#' \url{http://support.sas.com/kb/33/092.html}. The c is set to \eqn{c = 0.5}
#' if \eqn{r < \frac{2 n_0 n_1}{n_0 + n_1} + 1} and to \eqn{c = -0.5} if \eqn{r
#' > \frac{2 n_0 n_1}{n_0 + n_1} + 1}.
#' 
#' \bold{The Wald-Wolfowitz test} \verb{ } is a 2-sample nonparametric test to
#' evaluate if two continuous cumulative distributions are significantly
#' different or not. Ideally there should be no ties in the data. In practice
#' there is no problem with ties within a group, but if ties occur between
#' members of the different groups then there is no unique sequence of
#' observations. For example the data sets A: 10,14,17,19,34 and B:
#' 12,13,17,19,22 can give four possible sequences, with two possible values
#' for r (7 or 9). The "solution" to this is to list every possible
#' combination, and calculate the test statistic for each one. If all test
#' statistics are significant at the chosen level, then one can reject the null
#' hypothesis. If only some are significant, then Siegel (1956) suggests that
#' the average of the P-values is taken. Help for finding all permutations of
#' ties can be found at:
#' \url{https://stackoverflow.com/questions/47565066/all-possible-permutations-in-factor-variable-when-ties-exist-in-r}
#' 
#' However this solutions seems quite coarse and in general, the test should
#' not be used if there are more than one or two ties. We have better tests to
#' distinguish between two samples!
#' 
#' @name runsTest
#' @aliases runsTest runsTest.formula runsTest.default
#' @param x a dichotomous vector of data values or a (non-empty) numeric vector
#' of data values.
#' @param y an optional (non-empty) numeric vector of data values.
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
#' @param alternative a character string specifying the alternative hypothesis,
#' must be one of \code{"two.sided"} (default), \code{"less"} or
#' \code{"greater"}.
#' @param exact a logical indicating whether an exact p-value should be
#' computed. By default exact values will be calculated for small vectors with
#' a total length <= 30 and the normal approximation for longer ones.
#' @param correct a logical indicating whether to apply continuity correction
#' when computing the test statistic. Default is \code{TRUE}. Ignored if
#' \code{exact} is set to \code{TRUE}. See details.
#' @param na.rm defines if \code{NA}s should be omitted. Default is
#' \code{FALSE}.
#' @param \dots further arguments to be passed to or from methods.
#' @return A list with the following components.  \item{statistic}{z, the value
#' of the standardized runs statistic, if not exact p-values are computed.}
#' \item{parameter}{the number of runs, the total number of zeros (m) and ones
#' (n)} \item{p.value}{the p-value for the test.} \item{data.name}{a character
#' string giving the names of the data.} \item{alternative}{a character string
#' describing the alternative hypothesis.}
#' 
#' @seealso Run Length Encoding \code{\link{rle}}
#' @references Wackerly, D., Mendenhall, W. Scheaffer, R. L. (1986)
#' \emph{Mathematical Statistics with Applications}, 3rd Ed., Duxbury Press,
#' CA.
#' 
#' Wald, A. and Wolfowitz, J. (1940): On a test whether two samples are from
#' the same population, \emph{Ann. Math Statist}. 11, 147-162.
#' 
#' Siegel, S. (1956) \emph{Nonparametric Statistics for the Behavioural
#' Sciences}, McGraw-Hill Kogakusha, Tokyo.
#' 
#' @family topic.randomnessTests
#' @concept randomness
#' 
#' @examples
#' 
#' # x will be coerced to a dichotomous variable
#' x <- c("S","S", "T", "S", "T","T","T", "S", "T")
#' runsTest(x)
#' 
#' 
#' x <- c(13, 3, 14, 14, 1, 14, 3, 8, 14, 17, 9, 14, 13, 2, 16, 1, 3, 12, 13, 14)
#' runsTest(x)
#' # this will be treated as
#' runsTest(x > median(x))
#' 
#' plot( (x < median(x)) - 0.5, type="s", ylim=c(-1,1) )
#' abline(h=0)
#' 
#' set.seed(123)
#' x <- sample(0:1, size=100, replace=TRUE)
#' runsTest(x)
#' # As you would expect of values from a random number generator, the test fails to reject
#' # the null hypothesis that the data are random.
#' 
#' 
#' # SPSS example
#' x <- c(31,23,36,43,51,44,12,26,43,75,2,3,15,18,78,24,13,27,86,61,13,7,6,8)
#' runsTest(x, exact=TRUE)       # exact probability
#' runsTest(x, exact=FALSE)      # normal approximation
#' 
#' # SPSS example small dataset
#' x <- c(1, 1, 1, 1, 0, 0, 0, 0, 1, 1)
#' runsTest(x)
#' runsTest(x, exact=FALSE)
#' 
#' # if y is not NULL, the Wald-Wolfowitz-Test will be performed
#' A <- c(35,44,39,50,48,29,60,75,49,66)
#' B <- c(17,23,13,24,33,21,18,16,32)
#' 
#' runsTest(A, B, exact=TRUE)
#' runsTest(A, B, exact=FALSE)
#' 



#' @rdname runsTest
#' @export
runsTest <- function (x, ...)
  UseMethod("runsTest")



#' @rdname runsTest
#' @export
runsTest.formula <- local({
  
  # super elegant formula implementation
  # in fact we need nothing other, than is already implemented in
  # t.test.formula, besides the last call of zTest() instead of t.test()
  
  tf <- getS3method("t.test", "formula")
  
  new_body <- .replace_text_calls(body(tf), old="t.test", new="runsTest")
  
  new_fun <- tf
  body(new_fun) <- new_body
  
  new_fun
  
})



#' @rdname runsTest
#' @export
runsTest.default <- function(x, y=NULL, 
                             alternative=c("two.sided", "less", "greater"), 
                             exact=NULL, correct=TRUE, na.rm = FALSE, ...) {
 
  # exact values:
  # http://www.reiter1.com/Glossar/Wald_Wolfowitz.htm
  
  if(!is.null(y)) {
    dname <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
    
    xy <- cbind(c(x,y), c(rep(0, length(x)), rep(1, length(y))))
    xy <- xy[order(xy[,1], xy[,2]), 2]
    
    TIES <- length(intersect(x, y)) > 0
    
    res <- runsTest(x=xy, alternative=alternative, exact=exact, na.rm=na.rm)
    
    if (TIES)
      warning("cannot compute reliable p-values with inter-group ties")
    
    res$data.name <- dname
    res$method <- "Wald-Wolfowitz Runs Test"
    return(res)
  }
  
  alternative <- match.arg(alternative)
  dname <- deparse(substitute(x))
  
  if (na.rm) x <- na.omit(x)
  
  if(is.numeric(x) && length(unique(x)) > 2) {
    est <- median(x, na.rm=TRUE)
    names(est) <- "median(x)"
    x <- as.integer(x > est)
  } else {
    est <- NULL
  }
  
  x <- factor(x)
  if (nlevels(x) < 2)
    stop("Data must contain at least two distinct values")
  
  x <- as.numeric(x) - 1
  
  runs <- sum(diff(x) != 0) + 1
  m <- sum(x==0)
  n <- sum(x==1)
  
  if (is.null(exact))
    exact <- (m + n <= 30 && min(m,n) > 0)
  
  E <- 1 + 2*n*m / (n + m)
  s2 <- (2*n*m * (2*n*m - n - m)) / ((n + m)^2 * (n + m - 1))
  
  if (s2 <= 0) stop("Variance undefined")
  
    # this is the SPSS-Definition
    # http://publib.boulder.ibm.com/infocenter/spssstat/v20r0m0/index.jsp?topic=%2Fcom.ibm.spss.statistics.help%2Fidh_idd_npar_onesample_settings_tests_runs.htm

    statistic <- if(correct) {
    switch(as.character(cut(runs - E, c(-Inf,-0.5,0.5,Inf), 
                            labels=letters[1:3])),
           "a" = (runs - E + 0.5) / sqrt(s2),
           "b" = 0,
           "c" = (runs - E - 0.5) / sqrt(s2))
      
  } else {
    (runs - E) / sqrt(s2)
  }

  switch( alternative
          , "less" = {
            p.value <- if(exact) .pruns(runs, m, n,"less") else pnorm(statistic)
            alternative <- "true number of runs is less than expected"
          }
          , "greater" = {
            p.value = if(exact) .pruns(runs, m, n,"greater") else 1 - pnorm(statistic)
            alternative <- "true number of runs is greater than expected"
          }
          , "two.sided" = {
            p.value = if(exact) .pruns(runs,m,n,"two.sided") else 
                           2 * min(pnorm(statistic), 1 - pnorm(statistic))
            alternative <- "true number of runs is not equal the expected number"
          }
  )
  
  
  structure(list(
    statistic = if(!exact) c(z = statistic) else NULL,
    p.value = p.value,
    method = "Runs Test for Randomness",
    alternative = alternative,
    data.name = dname,
    # do not return estimate here ...
    parameter = c(runs=runs, m = m, n = n)
  ), class = "htest")
}





# == internal helper functions ==============================================


.pruns <- function(r, n1, n2, 
                   alternative=c("two.sided","less","greater")) {
  
  alternative <- match.arg(alternative)
  pruns_rcpp(r, n1, n2, alternative)
}


# replaced with cpp

# # function for calculating the denominator of the runs distribution
# .druns_nom <- function(r, n1, n2){
#   
#   pp <- vector(mode="numeric",length=length(r))
#   
#   for (i in seq_along(r)){
#     
#     if (2*r[i]%/%2==r[i]){
#       # even 2*k
#       k <- r[i]/2
#       pp[i] <- 2*choose(n1-1, k-1)*choose(n2-1, k-1)
#       
#     } else {
#       # odd 2*k+1
#       k <- (r[i]-1)/2
#       pp[i] <- choose(n1-1,k-1) * choose(n2-1,k) +
#         choose(n1-1,k)   * choose(n2-1,k-1)
#     }
#   }
#   
#   return(pp)
#   
# }
# 
# 
# 
# .pruns <- function(r, n1, n2, alternative=c("two.sided", "less", "greater")) {
#   
#   # source: randomizeBE
#   # author: D. Labes <detlewlabes at gmx.de>
#   
#   
#   alternative <- match.arg(alternative)
#   
#   n <- n1+n2
#   
#   if(r<=1) stop("Number of runs must be > 1")
#   if(r>n) stop("Number of runs must be < (n1+n2")
#   if(n1<1 | n2<1) return(0) #??? is not random!
#   
#   E <- 1 + 2*n1*n2/n
#   
#   denom <- choose(n,n1)
#   # how long should we make the r vector?
#   # in very unsymmetric cases only a few elements of
#   # pp = density have values > 0 if rmax=n1+n2
#   # number of runs possible: 2*m if n=m, 2*m+1 if m<n
#   rmax <- ifelse(n1==n2, 2*n1, 2*min(n1,n2)+1)
#   rv <- 2:rmax
#   pp <- .druns_nom(rv, n1, n2)
#   
#   # pL is p(R<=r) -> left/lower tail
#   pL <- sum(pp[rv<=r])/denom
#   #pU is p(R>=r) -> right/upper tail
#   pU <- 1 - sum(pp[rv<=(r-1)])/denom
#   
#   # Equn. 4.7 of the SPSS documentation
#   p2 <- sum(pp[abs(rv-E)>=abs(r-E)])/denom
#   
#   # Next is the rule from:
#   # Gibbons "Nonparametric Methods for Quantitative Analysis"
#   # 0.5 is to avoid p>1 if both pL and pU are >0.5
#   p2min <- 2*min(c(pL, pU, 0.5))
#   
#   # we are using the SPSS approach wich takes into account the
#   # unsymmetric form of the distribution if n1 << n2
#   
#   return(
#     switch( alternative
#             , "less" = pL
#             , "greater" = pU
#             , "two.sided" = p2
#     )
#   )
#   
# }
# 
