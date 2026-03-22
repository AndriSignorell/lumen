
#' The Triangular Distribution
#' 
#' The triangular distribution is a continuous distribution with a lower 
#' bound, an upper bound, and a mode, producing a piecewise linear, 
#' triangular-shaped density function. It is commonly used in risk 
#' assessment and simulation when only the minimum, maximum, and most 
#' likely value of a quantity are known.
#' 
#' Density, distribution function, quantile function, and random generation for
#' the triangular distribution with parameters \code{min}, \code{max}, and
#' \code{mode}.
#' 
#' Let \eqn{X} be a triangular random variable with parameters
#' \code{min=}\eqn{a}, \code{max=}\eqn{b}, and \code{mode=}\eqn{c}.
#' 
#' \emph{Probability Density and Cumulative Distribution Function} \cr The
#' density function of \eqn{X} is given by: \tabular{lll}{ \eqn{f(x; a, b, c)
#' =} \tab \eqn{\frac{2(x-a)}{(b-a)(c-a)}} \tab for \eqn{a \le x \le c} \cr
#' \tab \eqn{\frac{2(b-x)}{(b-a)(b-c)}} \tab for \eqn{c \le x \le b} \cr }
#' where \eqn{a < c < b}.
#' 
#' The cumulative distribution function of \eqn{X} is given by: \tabular{lll}{
#' \eqn{F(x; a, b, c) =} \tab \eqn{\frac{(x-a)^2}{(b-a)(c-a)}} \tab for \eqn{a
#' \le x \le c} \cr \tab \eqn{1 - \frac{(b-x)^2}{(b-a)(b-c)}} \tab for \eqn{c
#' \le x \le b} \cr } where \eqn{a < c < b}.
#' 
#' \emph{Quantiles} \cr The \eqn{p^th} quantile of \eqn{X} is given by:
#' \tabular{lll}{ \eqn{x_p =} \tab \eqn{a + \sqrt{(b-a)(c-a)p}} \tab for \eqn{0
#' \le p \le F(c)} \cr \tab \eqn{b - \sqrt{(b-a)(b-c)(1-p}} \tab for \eqn{F(c)
#' \le p \le 1} \cr } where \eqn{0 \le p \le 1}.
#' 
#' \emph{Random Numbers} \cr Random numbers are generated using the inverse
#' transformation method: \deqn{x = F^{-1}(u)} where \eqn{u} is a random
#' deviate from a uniform \eqn{[0, 1]} distribution.
#' 
#' \emph{Mean and Variance} \cr The mean and variance of \eqn{X} are given by:
#' \deqn{E(X) = \frac{a + b + c}{3}} \deqn{Var(X) = \frac{a^2 + b^2 + c^2 - ab
#' - ac - bc}{18}}
#' 
#' @name dtri
#' @aliases Triangular dtri ptri qTri rTri
#' @param x vector of quantiles.  Missing values (\code{NA}s) are allowed.
#' @param q vector of quantiles.  Missing values (\code{NA}s) are allowed.
#' @param p vector of probabilities between 0 and 1.  Missing values
#' (\code{NA}s) are allowed.
#' @param n sample size.  If \code{length(n)} is larger than 1, then
#' \code{length(n)} random values are returned.
#' @param min vector of minimum values of the distribution of the random
#' variable.  The default value is \code{min=0}.
#' @param max vector of maximum values of the random variable.  The default
#' value is \code{max=1}.
#' @param mode vector of modes of the random variable.  The default value is
#' \code{mode=1/2}.
#' @return \code{dtri} gives the density, \code{ptri} gives the distribution
#' function, \code{qTri} gives the quantile function, and \code{rTri} generates
#' random deviates.
#' @note The triangular distribution is so named because of the shape of its
#' probability density function.  The average of two independent identically
#' distributed uniform random variables with parameters \code{min=}\eqn{\alpha}
#' and \code{max=}\eqn{\beta} has a triangular distribution with parameters
#' \code{min=}\eqn{\alpha}, \code{max=}\eqn{\beta}, and
#' \code{mode=}\eqn{(\beta-\alpha)/2}.
#' 
#' The triangular distribution is sometimes used as an input distribution in
#' probability risk assessment.
#' 
#' @note
#' Based on code by Steven P. Millard. 
#' 
#' @seealso \link[stats:Uniform]{Uniform}, Probability Distributions and Random
#' Numbers.
#' 
#' @references Forbes, C., M. Evans, N. Hastings, and B. Peacock. (2011).
#' Statistical Distributions.  Fourth Edition. John Wiley and Sons, Hoboken,
#' NJ.
#' 
#' Johnson, N. L., S. Kotz, and N. Balakrishnan. (1995).  \emph{Continuous
#' Univariate Distributions, Volume 2}.  Second Edition. John Wiley and Sons,
#' New York.
#' 
#' @family topic.distributions
#' @concept continuous distribution
#' @concept triangular distribution
#' @concept bounded support
#' @concept dpqr
#' 
#' 
#' @examples
#' 
#'   # Density of a triangular distribution with parameters 
#'   # min=10, max=15, and mode=12, evaluated at 12, 13 and 14: 
#' 
#'   dtri(12:14, 10, 15, 12) 
#'   #[1] 0.4000000 0.2666667 0.1333333
#' 
#'   #----------
#' 
#'   # The cdf of a triangular distribution with parameters 
#'   # min=2, max=7, and mode=5, evaluated at 3, 4, and 5: 
#' 
#'   ptri(3:5, 2, 7, 5) 
#'   #[1] 0.06666667 0.26666667 0.60000000
#' 
#'   #----------
#' 
#'   # The 25'th percentile of a triangular distribution with parameters 
#'   # min=1, max=4, and mode=3: 
#' 
#'   qTri(0.25, 1, 4, 3) 
#'   #[1] 2.224745
#' 
#'   #----------
#' 
#'   # A random sample of 4 numbers from a triangular distribution with 
#'   # parameters min=3 , max=20, and mode=12. 
#'   # (Note: the call to set.seed simply allows you to reproduce this example.)
#' 
#'   set.seed(10) 
#'   rTri(4, 3, 20, 12) 
#'   #[1] 11.811593  9.850955 11.081885 13.539496
#' 

# Source: EnvStats
# author: Steven P. Millard (\email{EnvStats@ProbStatInfo.com})
# Version: 2.8.1



#' @rdname dtri
#' @export
dtri <- function (x, min = 0, max = 1, mode = 1/2) {
  names.x <- names(x)
  arg.mat <- .cbind.no.warn(x = as.vector(x), min = as.vector(min),
                           max = as.vector(max), mode = as.vector(mode))
  na.index <- .is_na_matrix(arg.mat)
  if (all(na.index))
    y <- rep(NA, nrow(arg.mat))
  else {
    y <- numeric(nrow(arg.mat))
    y[na.index] <- NA
    y.no.na <- y[!na.index]
    for (i in c("x", "min", "max", "mode")) assign(i, arg.mat[!na.index,
                                                              i])
    if (any(is.infinite(min)) || any(is.infinite(max)))
      stop("All non-missing values of 'min' and 'max' must be finite.")
    if (any(mode <= min) || any(max <= mode))
      stop(paste("All values of 'mode' must be larger than",
                 "the corresponding values of 'min', and all",
                 "values of 'max' must be larger than the", "corresponding values of 'mode'."))
    mmm <- max - min
    y.no.na <- 2 * ifelse(x <= mode, (x - min)/(mmm * (mode -
                                                         min)), (max - x)/(mmm * (max - mode)))
    y.no.na[y.no.na < 0] <- 0
    y[!na.index] <- y.no.na
  }
  if (!is.null(names.x))
    names(y) <- rep(names.x, length = length(y))
  else names(y) <- NULL
  y
}



#' @rdname dtri
#' @export
ptri <- function (q, min = 0, max = 1, mode = 1/2) {
  names.q <- names(q)
  arg.mat <- .cbind.no.warn(q = as.vector(q), min = as.vector(min),
                           max = as.vector(max), mode = as.vector(mode))
  na.index <- .is_na_matrix(arg.mat)
  if (all(na.index))
    p <- rep(NA, nrow(arg.mat))
  else {
    p <- numeric(nrow(arg.mat))
    p[na.index] <- NA
    p.no.na <- p[!na.index]
    for (i in c("q", "min", "max", "mode")) assign(i, arg.mat[!na.index,
                                                              i])
    if (any(is.infinite(min)) || any(is.infinite(max)))
      stop("All non-missing values of 'min' and 'max' must be finite.")
    if (any(mode <= min) || any(max <= mode))
      stop(paste("All values of 'mode' must be larger than",
                 "the corresponding values of 'min', and all",
                 "values of 'max' must be larger than the", "corresponding values of 'mode'."))
    q.low <- q <= min
    p.no.na[q.low] <- 0
    q.high <- q >= max
    p.no.na[q.high] <- 1
    if (any(index <- !(q.low | q.high))) {
      for (i in c("q", "min", "max", "mode")) assign(i,
                                                     get(i)[index])
      mmm <- max - min
      p.no.na[index] <- ifelse(q <= mode, (q - min)^2/(mmm *
                                                         (mode - min)), 1 - ((max - q)^2/(mmm * (max -
                                                                                                   mode))))
    }
    p[!na.index] <- p.no.na
  }
  if (!is.null(names.q))
    names(p) <- rep(names.q, length = length(p))
  else names(p) <- NULL
  p
}



#' @rdname dtri
#' @export
qTri <- function (p, min = 0, max = 1, mode = 1/2) {
  names.p <- names(p)
  arg.mat <- .cbind.no.warn(p = as.vector(p), min = as.vector(min),
                           max = as.vector(max), mode = as.vector(mode))
  na.index <- .is_na_matrix(arg.mat)
  if (all(na.index))
    q <- rep(NA, nrow(arg.mat))
  else {
    q <- numeric(nrow(arg.mat))
    q[na.index] <- NA
    q.no.na <- q[!na.index]
    for (i in c("p", "min", "max", "mode")) assign(i, arg.mat[!na.index,
                                                              i])
    if (any(p < 0) || any(p > 1))
      stop("All non-missing values of 'p' must be between 0 and 1.")
    if (any(is.infinite(min)) || any(is.infinite(max)))
      stop("All non-missing values of 'min' and 'max' must be finite.")
    if (any(mode <= min) || any(max <= mode))
      stop(paste("All values of 'mode' must be larger than",
                 "the corresponding values of 'min', and all",
                 "values of 'max' must be larger than the", "corresponding values of 'mode'."))
    q.no.na[p == 0] <- min[p == 0]
    q.no.na[p == 1] <- max[p == 1]
    if (any(index <- 0 < p & p < 1)) {
      for (i in c("p", "min", "max", "mode")) assign(i,
                                                     get(i)[index])
      mmm <- max - min
      q.no.na[index] <- ifelse(p <= ptri(mode, min = min,
                                         max = max, mode = mode), min + sqrt(mmm * (mode -
                                                                                      min) * p), max - sqrt(mmm * (max - mode) * (1 -
                                                                                                                                    p)))
    }
    q[!na.index] <- q.no.na
  }
  if (!is.null(names.p))
    names(q) <- rep(names.p, length = length(q))
  else names(q) <- NULL
  q
}


#' @rdname dtri
#' @export
rTri <- function (n, min = 0, max = 1, mode = 1/2) {
  ln <- length(n)
  if (ln < 1)
    stop("'n' must be non-empty.")
  if (ln > 1)
    n <- ln
  else {
    if (is.na(n) || n <= 0 || n != trunc(n))
      stop("'n' must be a positive integer or vector.")
  }
  arg.mat <- .cbind.no.warn(dum = rep(1, n), min = as.vector(min),
                           max = as.vector(max), mode = as.vector(mode))[, -1, drop = FALSE]
  if (n < nrow(arg.mat))
    arg.mat <- arg.mat[1:n, , drop = FALSE]
  for (i in c("min", "max", "mode")) assign(i, arg.mat[, i])
  na.index <- .is_na_matrix(arg.mat)
  if (all(na.index))
    return(rep(NA, n))
  else {
    if (any(is.infinite(min)) || any(is.infinite(max)))
      stop("All non-missing values of 'min' and 'max' must be finite.")
    if (any(mode <= min) || any(max <= mode))
      stop(paste("All values of 'mode' must be larger than",
                 "the corresponding values of 'min', and all",
                 "values of 'max' must be larger than the", "corresponding values of 'mode'."))
    return(qTri(p = runif(n), min = min, max = max, mode = mode))
  }
}


# == internal helper functions ======================================================


.is_na_matrix <- function (mat, rows = TRUE) {
  
  if (!is.matrix(mat))
    stop("'mat' must be a matrix or data frame.")
  if (rows)
    return(apply(mat, 1, function(x) any(is.na(x))))
  else return(apply(mat, 2, function(x) any(is.na(x))))
}


.cbind.no.warn <- function (..., deparse.level = 1) {
  
  oldopts <- options(warn = -1)
  on.exit(options(oldopts))
  base::cbind(..., deparse.level = deparse.level)
}



