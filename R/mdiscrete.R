#' Mean and Variance of Discrete Distributions
#'
#' Computes the mean and variance of common discrete distributions given 
#' their parameters.
#'
#' @details
#' \tabular{lll}{
#'   \strong{Distribution\verb{       } } \tab \strong{Mean\verb{                 } } \tab **Variance** \cr
#'   Binomial \tab
#'     \eqn{np} \tab
#'     \eqn{np(1-p)} \cr
#'   Poisson \tab
#'     \eqn{\lambda} \tab
#'     \eqn{\lambda} \cr
#'   Geometric \tab
#'     \eqn{\frac{1-p}{p}} \tab
#'     \eqn{\frac{1-p}{p^2}} \cr
#'   Negative binomial \tab
#'     \eqn{\frac{r(1-p)}{p}} \tab
#'     \eqn{\frac{r(1-p)}{p^2}} \cr
#'   Hypergeometric \tab
#'     \eqn{\frac{km}{N}} \tab
#'     \eqn{\frac{km}{N}\frac{n}{N}\frac{N-k}{N-1}} \cr
#'   Benford \tab
#'     \eqn{\sum_d d\log_{10}\left(1 + \frac{1}{d}\right)} \tab
#'     \eqn{\sum_d d^2\log_{10}\left(1 + \frac{1}{d}\right) - \mu^2} \cr
#' }
#'
#' For the binomial distribution, \eqn{n} = `size`; for the negative
#' binomial distribution, \eqn{r} = `size`; and for the hypergeometric
#' distribution, \eqn{N = m + n}. For Benford's distribution, the sum runs
#' over \eqn{d \in \{1,\ldots,9\}} for `nDigits = 1` and
#' \eqn{d \in \{10,\ldots,99\}} for `nDigits = 2`. As there is no
#' closed-form solution, the moments are computed numerically.
#' 
#' @param size number of trials (binomial, negative binomial).
#' @param prob probability of success on each trial (binomial, geometric, 
#'   negative binomial).
#' @param lambda mean (Poisson).
#' @param m number of white balls in the urn (hypergeometric).
#' @param n number of black balls in the urn (hypergeometric).
#' @param k number of balls drawn (hypergeometric).
#' @param nDigits number of leading digits for Benford's distribution,
#'   either `1` (default, support \{1,...,9\}) or `2`
#'   (support \{10,...,99\}).
#'
#' @return A named numeric vector with elements `mean` and 
#'   `variance`.
#'
#' @references
#' Forbes, C., Evans, M., Hastings, N. and Peacock, B. (2011)
#' *Statistical Distributions*. Fourth Edition. Wiley.
#'
#' Johnson, N. L., Kotz, S. and Balakrishnan, N. (1995)
#' *Continuous Univariate Distributions*, Vol. 2. Wiley.
#' 
#' @seealso [Binomial()], [Poisson()],
#'   [Geometric()], [NegBinomial()],
#'   [Hypergeometric()], [distributions-overview]
#'
#' @concept distribution-summary
#' @concept moment
#'
#' @examples
#' mbinom(size = 10, prob = 0.5)
#' mpois(lambda = 3)
#' mgeom(prob = 0.3)
#' mnbinom(size = 5, prob = 0.3)
#' mhyper(m = 10, n = 5, k = 4)
#' mbenford(nDigits = 1)
#' mbenford(nDigits = 2)
#'
#' @name disc.moments
NULL

#' @rdname disc.moments
#' @export
mbinom <- function(size, prob) {

  .assertScalar(size, lower = 0, integerValued = TRUE)
  .assertScalar(prob, lower = 0, upper = 1)

  c(mean     = size * prob,
    variance = size * prob * (1 - prob))
}

#' @rdname disc.moments
#' @export
mpois <- function(lambda) {

  .assertScalar(lambda, lower = 0)

  c(mean     = lambda,
    variance = lambda)
}

#' @rdname disc.moments
#' @export
mgeom <- function(prob) {

  .assertScalar(prob, lower = 0, upper = 1, strictLower = TRUE)

  c(mean     = (1 - prob) / prob,
    variance = (1 - prob) / prob^2)
}

#' @rdname disc.moments
#' @export
mnbinom <- function(size, prob) {

  # 'size' need not be a whole number: dnbinom() admits the continuous
  # gamma-mixture parametrisation as well
  .assertScalar(size, lower = 0, strictLower = TRUE)
  .assertScalar(prob, lower = 0, upper = 1, strictLower = TRUE)

  c(mean     = size * (1 - prob) / prob,
    variance = size * (1 - prob) / prob^2)
}

#' @rdname disc.moments
#' @export
mhyper <- function(m, n, k) {

  .assertScalar(m, lower = 0, integerValued = TRUE)
  .assertScalar(n, lower = 0, integerValued = TRUE)

  N <- m + n

  # the variance divides by N - 1
  if (N < 2)
    stop("'m + n' must be at least 2", call. = FALSE)

  .assertScalar(k, lower = 0, upper = N, integerValued = TRUE)

  c(mean     = k * m / N,
    variance = k * m / N * n / N * (N - k) / (N - 1))
}

#' @rdname disc.moments
#' @export
mbenford <- function(nDigits = 1) {

  .assertScalar(nDigits, lower = 1, upper = 2, integerValued = TRUE)

  d <- if (nDigits == 1) 1:9 else 10:99
  p <- log10(1 + 1/d)
  mu <- sum(d * p)
  c(mean     = mu,
    variance = sum(d^2 * p) - mu^2)
}



