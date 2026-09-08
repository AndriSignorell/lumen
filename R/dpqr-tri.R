
#' Triangular Distribution
#' 
#' The triangular distribution is a continuous distribution with a lower 
#' bound, an upper bound, and a mode, producing a piecewise linear, 
#' triangular-shaped density function. It is commonly used in risk 
#' assessment and simulation when only the minimum, maximum, and most 
#' likely value of a quantity are known.
#' 
#' Density, distribution function, quantile function, and random generation for
#' the triangular distribution with parameters `min`, `max`, and
#' `mode`.
#' 
#' Let \eqn{X} be a triangular random variable with parameters
#' `min=`\eqn{a}, `max=`\eqn{b}, and `mode=`\eqn{c}.
#' 
#' *Probability Density and Cumulative Distribution Function* \cr The
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
#' *Quantiles* \cr The \eqn{p^th} quantile of \eqn{X} is given by:
#' \tabular{lll}{ \eqn{x_p =} \tab \eqn{a + \sqrt{(b-a)(c-a)p}} \tab for \eqn{0
#' \le p \le F(c)} \cr \tab \eqn{b - \sqrt{(b-a)(b-c)(1-p}} \tab for \eqn{F(c)
#' \le p \le 1} \cr } where \eqn{0 \le p \le 1}.
#' 
#' *Random Numbers* \cr Random numbers are generated using the inverse
#' transformation method: \deqn{x = F^{-1}(u)} where \eqn{u} is a random
#' deviate from a uniform \eqn{[0, 1]} distribution.
#' 
#' *Mean and Variance* \cr The mean and variance of \eqn{X} are given by:
#' \deqn{E(X) = \frac{a + b + c}{3}} \deqn{Var(X) = \frac{a^2 + b^2 + c^2 - ab
#' - ac - bc}{18}}
#' 
#' The triangular distribution is so named because of the shape of its
#' probability density function. The average of two independent identically
#' distributed uniform random variables with parameters `min=`\eqn{\alpha}
#' and `max=`\eqn{\beta} has a triangular distribution with parameters
#' `min=`\eqn{\alpha}, `max=`\eqn{\beta}, and
#' `mode=`\eqn{(\alpha+\beta)/2}.
#' 
#' @name dpqr-tri
#' @aliases Triangular dtri ptri qtri rtri
#' 
#' @param x vector of quantiles.  Missing values (`NA`s) are allowed.
#' @param q vector of quantiles.  Missing values (`NA`s) are allowed.
#' @param p vector of probabilities between 0 and 1.  Missing values
#' (`NA`s) are allowed.
#' @param n sample size.  If `length(n)` is larger than 1, then
#' `length(n)` random values are returned.
#' @param min vector of minimum values of the distribution of the random
#' variable.  The default value is `min=0`.
#' @param max vector of maximum values of the random variable.  The default
#' value is `max=1`.
#' @param mode vector of modes of the random variable.  The default value is
#' `mode=1/2`.  The parameters must satisfy \eqn{min < mode < max}.
#' @param log,log.p logical; if `TRUE`, probabilities `p` are given as
#' `log(p)` and the density is returned on the log scale.
#' @param lower.tail logical; if `TRUE` (default), probabilities are 
#' \verb{P[X <= x]}, otherwise, P\verb{[X > x]}.
#' @return `dtri()` gives the density, `ptri()` gives the
#' distribution function, `qtri()` gives the quantile function, and
#' `rtri()` generates random deviates.
#' 
#' @details
#' The triangular distribution is sometimes used as an input distribution in
#' probability risk assessment.
#' 
#' @note
#' Based on code by Steven P. Millard previously published in
#' the \pkg{EnvStats} package, adapted to conform to package standards.
#' 
#' @seealso [distributions-overview]; [Uniform][Uniform]
#' 
#' @references Forbes, C., M. Evans, N. Hastings, and B. Peacock. (2011).
#' Statistical Distributions.  Fourth Edition. John Wiley and Sons, Hoboken,
#' NJ.
#' 
#' Johnson, N. L., S. Kotz, and N. Balakrishnan. (1995).  *Continuous
#' Univariate Distributions, Volume 2*.  Second Edition. John Wiley and Sons,
#' New York.
#' 
#' 
#' @examples
#' 
#' # Density of a triangular distribution with parameters 
#' # min=10, max=15, and mode=12, evaluated at 12, 13 and 14: 
#' 
#' dtri(12:14, 10, 15, 12) 
#' ## [1] 0.4000000 0.2666667 0.1333333
#' 
#' # The cdf of a triangular distribution with parameters 
#' # min=2, max=7, and mode=5, evaluated at 3, 4, and 5: 
#' 
#' ptri(3:5, 2, 7, 5) 
#' ## [1] 0.06666667 0.26666667 0.60000000
#' 
#' # The 25'th percentile of a triangular distribution with parameters 
#' # min=1, max=4, and mode=3: 
#' 
#' qtri(0.25, 1, 4, 3) 
#' ## [1] 2.224745
#' 
#' # A random sample of 4 numbers from a triangular distribution with 
#' # parameters min=3 , max=20, and mode=12. 
#' # (Note: the call to set.seed simply allows you to reproduce this example.)
#' 
#' set.seed(10) 
#' rtri(4, 3, 20, 12) 
#' ## [1] 11.811593  9.850955 11.081885 13.539496
#' 

# Source: EnvStats
# author: Steven P. Millard (\email{EnvStats@ProbStatInfo.com})
# Version: 2.8.1



#' @rdname dpqr-tri
#' @concept distribution-function
#' @concept sampling
#' @export
dtri <- function (x, min = 0, max = 1, mode = 1/2, log = FALSE) {

  # 'min' and 'max' shadow the base functions of the same name
  a <- min; b <- max; m <- mode
  .checkTri(a, b, m)

  nms <- names(x)
  n <- base::max(length(x), length(a), length(b), length(m))
  x <- rep_len(x, n); a <- rep_len(a, n)
  b <- rep_len(b, n); m <- rep_len(m, n)

  d   <- numeric(n)
  bad <- is.na(x) | is.na(a) | is.na(b) | is.na(m)

  lo <- which(!bad & x >= a & x <= m)
  hi <- which(!bad & x >  m & x <= b)
  d[lo] <- 2 * (x[lo] - a[lo]) / ((b[lo] - a[lo]) * (m[lo] - a[lo]))
  d[hi] <- 2 * (b[hi] - x[hi]) / ((b[hi] - a[hi]) * (b[hi] - m[hi]))
  d[bad] <- NA

  if (!is.null(nms)) names(d) <- rep_len(nms, n)
  if (log) base::log(d) else d
}



#' @rdname dpqr-tri
#' @export
ptri <- function (q, min = 0, max = 1, mode = 1/2, lower.tail = TRUE,
                  log.p = FALSE) {

  a <- min; b <- max; m <- mode
  .checkTri(a, b, m)

  nms <- names(q)
  n <- base::max(length(q), length(a), length(b), length(m))
  q <- rep_len(q, n); a <- rep_len(a, n)
  b <- rep_len(b, n); m <- rep_len(m, n)

  p   <- numeric(n)
  bad <- is.na(q) | is.na(a) | is.na(b) | is.na(m)

  lo <- which(!bad & q >  a & q <= m)
  hi <- which(!bad & q >  m & q <  b)
  p[which(!bad & q >= b)] <- 1
  p[lo] <- (q[lo] - a[lo])^2 / ((b[lo] - a[lo]) * (m[lo] - a[lo]))
  p[hi] <- 1 - (b[hi] - q[hi])^2 / ((b[hi] - a[hi]) * (b[hi] - m[hi]))
  p[bad] <- NA

  if (!lower.tail) p <- 1 - p
  if (!is.null(nms)) names(p) <- rep_len(nms, n)
  if (log.p) base::log(p) else p
}



#' @rdname dpqr-tri
#' @export
qtri <- function (p, min = 0, max = 1, mode = 1/2, lower.tail = TRUE,
                  log.p = FALSE) {

  a <- min; b <- max; m <- mode
  .checkTri(a, b, m)

  nms <- names(p)
  p <- .qProb(p, lower.tail = lower.tail, log.p = log.p)
  n <- base::max(length(p), length(a), length(b), length(m))
  p <- rep_len(p, n); a <- rep_len(a, n)
  b <- rep_len(b, n); m <- rep_len(m, n)

  # F(mode), the probability at which the two branches meet
  pm <- (m - a) / (b - a)
  q  <- ifelse(p <= pm,
               a + sqrt((b - a) * (m - a) * p),
               b - sqrt((b - a) * (b - m) * (1 - p)))

  # ifelse() turns a NaN condition into NA; keep the two apart
  q[is.na(p)] <- p[is.na(p)]

  if (!is.null(nms)) names(q) <- rep_len(nms, n)
  q
}


#' @rdname dpqr-tri
#' @export
rtri <- function (n, min = 0, max = 1, mode = 1/2) {

  if (length(n) > 1) n <- length(n)
  .assertScalar(n, lower = 1, integerValued = TRUE)

  qtri(runif(n), min = min, max = max, mode = mode)
}


# == internal helper functions ===============================================


# the three parameters must satisfy min < mode < max; missing values are
# passed over here and propagate to the result
#
# @noRd
.checkTri <- function(min, max, mode) {

  n <- base::max(length(min), length(max), length(mode))
  a <- rep_len(min, n); b <- rep_len(max, n); m <- rep_len(mode, n)
  ok <- !(is.na(a) | is.na(b) | is.na(m))

  if (!is.numeric(a) || !is.numeric(b) || !is.numeric(m))
    stop("'min', 'max' and 'mode' must be numeric", call. = FALSE)

  if (any(!is.finite(a[ok])) || any(!is.finite(b[ok])))
    stop("'min' and 'max' must be finite", call. = FALSE)

  if (any(a[ok] >= m[ok]) || any(m[ok] >= b[ok]))
    stop("'min', 'mode' and 'max' must satisfy min < mode < max",
         call. = FALSE)

  invisible(TRUE)
}
