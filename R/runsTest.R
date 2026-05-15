
#' Runs Test for Randomness
#'
#' A nonparametric test for randomness of a sequence, based on the number
#' of runs (consecutive sequences of identical values or values above/below
#' a threshold).
#'
#' Performs a test whether the elements of \code{x} are serially independent
#' by counting how many runs there are above and below a threshold. If
#' \code{y} is supplied, a two-sample Wald-Wolfowitz test for the equality
#' of two distributions is computed.
#'
#' \bold{The runs test for randomness} requires a dichotomous sequence. For
#' a numeric variable \code{x} with more than two distinct values, the
#' sequence is dichotomised by comparing each observation to the median.
#' Observations exactly equal to the median are removed before the test, as
#' is standard in the runs test literature. To use a different threshold,
#' pass a logical vector directly: \code{runsTest(x > mean(x))}.
#'
#' The normal approximation uses the expected number of runs under the null
#' \deqn{\mu_r = \frac{2 n_0 n_1}{n_0 + n_1} + 1}
#' and its variance
#' \deqn{\sigma_r^2 = \frac{2 n_0 n_1 (2 n_0 n_1 - n_0 - n_1)}{
#'   (n_0 + n_1)^2 (n_0 + n_1 - 1)}}
#' as \eqn{\hat{z} = (r - \mu_r + c) / \sigma_r}, where \eqn{n_0, n_1} are
#' the number of values below/above the threshold and \eqn{r} is the number
#' of runs.
#'
#' Setting \code{correct = TRUE} applies a continuity correction as SAS (and
#' SPSS for \eqn{n < 50}) does: \eqn{c = 0.5} if \eqn{r < \mu_r} and
#' \eqn{c = -0.5} if \eqn{r > \mu_r}.
#'
#' The exact p-value is computed from the conditional distribution of the
#' number of runs given \eqn{n_0} and \eqn{n_1}, implemented in C++.
#'
#' \bold{Interpretation of alternatives:}
#' \itemize{
#'   \item \code{"less"}: fewer runs than expected — indicates clustering
#'     or positive serial correlation.
#'   \item \code{"greater"}: more runs than expected — indicates alternation
#'     or negative serial correlation.
#'   \item \code{"two.sided"}: departure in either direction.
#' }
#'
#' \bold{The Wald-Wolfowitz test} is a two-sample nonparametric test for the
#' equality of two continuous distributions against general alternatives.
#' Exact p-values are not valid in the presence of inter-group ties; a
#' warning is issued when such ties are detected.
#'
#' @name runsTest
#' @aliases runsTest runsTest.formula runsTest.default
#'
#' @param x a dichotomous vector of data values or a (non-empty) numeric
#'   vector of data values.
#' @param y an optional (non-empty) numeric vector of data values for the
#'   Wald-Wolfowitz two-sample test.
#' @param formula a formula of the form \code{lhs ~ rhs} where \code{lhs}
#'   gives the data values and \code{rhs} the corresponding groups (exactly
#'   two groups required).
#' @param data an optional data frame containing the variables in
#'   \code{formula}.
#' @param subset an optional vector specifying a subset of observations to
#'   be used.
#' @param na.action a function which indicates what should happen when the
#'   data contain \code{NA}s. Defaults to \code{getOption("na.action")}.
#' @param alternative a character string specifying the alternative
#'   hypothesis, must be one of \code{"two.sided"} (default),
#'   \code{"less"} (fewer runs, clustering) or \code{"greater"} (more
#'   runs, alternation).
#' @param exact a logical indicating whether an exact p-value should be
#'   computed. By default exact values are calculated for
#'   \eqn{n_0 + n_1 \le 30} and the normal approximation otherwise.
#' @param correct a logical indicating whether to apply a continuity
#'   correction when computing the test statistic. Default is \code{TRUE}.
#'   Ignored when \code{exact = TRUE}.
#' @param na.rm logical. If \code{TRUE} (default), \code{NA}s are removed
#'   before the test.
#' @param \dots further arguments passed to methods.
#'
#' @return A list of class \code{"htest"} containing:
#' \item{statistic}{the standardized z-statistic (normal approximation only;
#'   \code{NULL} for exact tests).}
#' \item{parameter}{named vector with the number of runs, \code{m}
#'   (\eqn{n_0}: observations below threshold) and \code{n}
#'   (\eqn{n_1}: observations above threshold).}
#' \item{p.value}{the p-value for the test.}
#' \item{alternative}{a character string describing the alternative
#'   hypothesis.}
#' \item{method}{a character string indicating the test performed.}
#' \item{data.name}{a character string giving the name of the data.}
#'
#' @seealso \code{\link{rle}}, \code{\link{bartelsRankTest}},
#'   \code{\link{vonNeumannTest}}
#'
#' @references
#' Wackerly, D., Mendenhall, W., Scheaffer, R. L. (1986)
#' \emph{Mathematical Statistics with Applications}, 3rd Ed.,
#' Duxbury Press, CA.
#'
#' Wald, A. and Wolfowitz, J. (1940): On a test whether two samples are
#' from the same population. \emph{Annals of Mathematical Statistics}
#' \bold{11}, 147--162.
#'
#' Siegel, S. (1956) \emph{Nonparametric Statistics for the Behavioural
#' Sciences}, McGraw-Hill Kogakusha, Tokyo.
#'
#' @examples
#' # dichotomous character vector
#' x <- c("S","S","T","S","T","T","T","S","T")
#' runsTest(x)
#'
#' # numeric vector: compared against median; ties at median removed
#' x <- c(13,3,14,14,1,14,3,8,14,17,9,14,13,2,16,1,3,12,13,14)
#' runsTest(x)
#'
#' # SPSS example
#' x <- c(31,23,36,43,51,44,12,26,43,75,2,3,15,18,78,24,13,27,86,61,13,7,6,8)
#' runsTest(x, exact = TRUE)
#' runsTest(x, exact = FALSE)
#'
#' # Wald-Wolfowitz two-sample test
#' A <- c(35,44,39,50,48,29,60,75,49,66)
#' B <- c(17,23,13,24,33,21,18,16,32)
#' runsTest(A, B, exact = TRUE)
#'
#' @rdname runsTest
#' @family test.correlation
#' @concept hypothesis-testing
#' @concept nonparametric
#' @concept time-series


#' @export
runsTest <- function(x, ...) UseMethod("runsTest")



#' @rdname runsTest
#' @export
runsTest.formula <- function(formula,
                             data,
                             subset,
                             na.action = na.pass,
                             ...) {
  
  if (missing(formula) || length(formula) != 3L)
    stop("'formula' missing or incorrect")
  
  args <- list(
    formula   = formula,
    na.action = na.action,
    allowed   = "two.sample.independent"
  )
  
  if (!missing(data))
    args$data <- data
  
  if (!missing(subset))
    args$subset <- substitute(subset)
  
  d <- do.call(bedrock::resolveFormula, args)
  
  res <- runsTest.default(
    x = d$x,
    y = d$y,
    ...
  )
  
  res$data.name <- d$data.name
  res
}


#' @rdname runsTest
#' @export
runsTest.default <- function(x,
                             y           = NULL,
                             alternative = c("two.sided", "less", "greater"),
                             exact       = NULL,
                             correct     = TRUE,
                             na.rm       = TRUE,
                             ...) {
  
  ## -------------------------------------------------------------------
  ## Two-sample Wald-Wolfowitz
  ## -------------------------------------------------------------------
  if (!is.null(y)) {
    
    if (length(x) == 0L || length(y) == 0L)
      stop("'x' and 'y' must be non-empty")
    
    DNAME <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
    
    grp <- c(rep(0L, length(x)), rep(1L, length(y)))
    val <- c(x, y)
    
    ## Inter-group ties: any value appearing in both groups
    if (any(val[grp == 0L] %in% val[grp == 1L]))
      warning(
        "inter-group ties detected; ",
        "exact p-values are not valid in this case"
      )
    
    ord <- order(val, grp)
    xy  <- grp[ord]
    
    res           <- runsTest.default(
      x           = xy,
      alternative = alternative,
      exact       = exact,
      na.rm       = FALSE
    )
    res$data.name <- DNAME
    res$method    <- "Wald-Wolfowitz Runs Test"
    return(res)
  }
  
  ## -------------------------------------------------------------------
  ## One-sample runs test
  ## -------------------------------------------------------------------
  DNAME       <- deparse(substitute(x))
  alternative <- match.arg(alternative)
  
  if (na.rm) {
    x <- x[!is.na(x)]
  } else if (anyNA(x)) {
    stop("'x' contains NA values; set na.rm = TRUE to remove them")
  }
  
  if (length(x) == 0L)
    stop("'x' is empty after removing NA values")
  
  ## Dichotomise numeric input: compare to median, remove median ties
  if (is.numeric(x) && length(unique(x)) > 2L) {
    s <- sign(x - median(x))
    x <- as.integer(s[s != 0L] > 0L)
    
  } else {
    ux <- unique(x)
    if (length(ux) < 2L)
      stop("data must contain at least two distinct values")
    if (length(ux) > 2L)
      stop("data must contain at most two distinct values")
    x <- as.integer(factor(x)) - 1L
  }
  
  
  if (length(unique(x)) < 2L)
    stop(
      "all observations fall on one side of the threshold; ",
      "the test is not defined"
    )
  
  runs <- sum(diff(x) != 0L) + 1L
  n0   <- sum(x == 0L)
  n1   <- sum(x == 1L)
  
  if (is.null(exact))
    exact <- (n0 + n1 <= 30L && min(n0, n1) > 0L)
  
  mu     <- 1 + 2 * n0 * n1 / (n0 + n1)
  sigma2 <- (2 * n0 * n1 * (2 * n0 * n1 - n0 - n1)) /
    ((n0 + n1)^2 * (n0 + n1 - 1L))
  
  if (sigma2 <= 0)
    stop("variance is undefined (all observations in one group?)")
  
  ## Continuity correction (SPSS/SAS convention)
  cc <- if (correct && !exact) {
    if      (runs - mu < -0.5)  0.5
    else if (runs - mu >  0.5) -0.5
    else                         0
  } else 0
  
  z <- (runs - mu + cc) / sqrt(sigma2)
  
  PVAL <- if (exact) {
    .pruns(runs, n0, n1, alternative)
  } else {
    switch(
      alternative,
      "less"      = pnorm(z),
      "greater"   = pnorm(z, lower.tail = FALSE),
      "two.sided" = 2 * min(pnorm(z), pnorm(z, lower.tail = FALSE))
    )
  }
  
  ALT <- switch(
    alternative,
    "less"      = "true number of runs is less than expected (clustering)",
    "greater"   = "true number of runs is greater than expected (alternation)",
    "two.sided" = "true number of runs is not equal to the expected number"
  )
  
  structure(
    list(
      statistic   = if (!exact) c(z = z) else NULL,
      parameter   = c(runs = runs, m = n0, n = n1),
      p.value     = as.numeric(PVAL),
      alternative = ALT,
      method      = "Runs Test for Randomness",
      data.name   = DNAME
    ),
    class = "htest"
  )
}



# == internal helper functions ================================================


# exact p-value via C++
.pruns <- function(r, n1, n2,
                   alternative = c("two.sided", "less", "greater")) {
  alternative <- match.arg(alternative)
  pruns_rcpp(r, n1, n2, alternative)
}





