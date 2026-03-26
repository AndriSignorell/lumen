#' Mean and Variance of Discrete Distributions
#'
#' Computes the mean and variance of common discrete distributions given 
#' their parameters.
#'
#' @details
#' \strong{Binomial:}\cr
#' \eqn{\mu = n \cdot p}\cr
#' \eqn{\mathrm{Var}(X) = n \cdot p \cdot (1-p)}
#'
#' \strong{Poisson:}\cr
#' \eqn{\mu = \lambda}\cr
#' \eqn{\mathrm{Var}(X) = \lambda}
#'
#' \strong{Geometric:}\cr
#' \eqn{\mu = \frac{1-p}{p}}\cr
#' \eqn{\mathrm{Var}(X) = \frac{1-p}{p^2}}
#'
#' \strong{Negative Binomial:}\cr
#' \eqn{\mu = \frac{r(1-p)}{p}}\cr
#' \eqn{\mathrm{Var}(X) = \frac{r(1-p)}{p^2}}
#'
#' \strong{Hypergeometric:}\cr
#' \eqn{\mu = \frac{km}{N}}\cr
#' \eqn{\mathrm{Var}(X) = \frac{km}{N} \cdot \frac{n}{N} \cdot \frac{N-k}{N-1}}
#' 
#' \strong{Benford:}\cr
#' \eqn{\mu = \sum_{d} d \cdot \log_{10}\left(1 + \frac{1}{d}\right)}\cr
#' \eqn{\mathrm{Var}(X) = \sum_{d} d^2 \cdot \log_{10}\left(1 + 
#'   \frac{1}{d}\right) - \mu^2}\cr
#' where the sum runs over \eqn{d \in \{1,\ldots,9\}} for \code{ndigits = 1}
#' and \eqn{d \in \{10,\ldots,99\}} for \code{ndigits = 2}. As there is no
#' closed-form solution, the moments are computed numerically.
#' 
#' @param size number of trials (binomial, negative binomial).
#' @param prob probability of success on each trial (binomial, geometric, 
#'   negative binomial).
#' @param lambda mean (Poisson).
#' @param m number of white balls in the urn (hypergeometric).
#' @param n number of black balls in the urn (hypergeometric).
#' @param k number of balls drawn (hypergeometric).
#' @param ndigits number of leading digits for Benford's distribution,
#'   either \code{1} (default, support \{1,...,9\}) or \code{2}
#'   (support \{10,...,99\}).
#'
#' @return A named numeric vector with elements \code{mean} and 
#'   \code{variance}.
#'
#' @references
#' Forbes, C., Evans, M., Hastings, N. and Peacock, B. (2011)
#' \emph{Statistical Distributions}. Fourth Edition. Wiley.
#'
#' Johnson, N. L., Kotz, S. and Balakrishnan, N. (1995)
#' \emph{Continuous Univariate Distributions}, Vol. 2. Wiley.
#' 
#' @seealso \code{\link[stats]{Binomial}}, \code{\link[stats]{Poisson}},
#'   \code{\link[stats]{Geometric}}, \code{\link[stats]{NegBinomial}},
#'   \code{\link[stats]{Hypergeometric}}
#'
#' @family topic.distributions.moments
#' @concept discrete distribution
#' @concept moments
#'
#' @examples
#' mbinom(size = 10, prob = 0.5)
#' mpois(lambda = 3)
#' mgeom(prob = 0.3)
#' mnbinom(size = 5, prob = 0.3)
#' mhyper(m = 10, n = 5, k = 4)
#' mbenford(ndigits = 1)
#' mbenford(ndigits = 2)
#'
#' @name disc.moments
NULL

#' @rdname disc.moments
#' @export
mbinom <- function(size, prob) {
  c(mean     = size * prob,
    variance = size * prob * (1 - prob))
}

#' @rdname disc.moments
#' @export
mpois <- function(lambda) {
  c(mean     = lambda,
    variance = lambda)
}

#' @rdname disc.moments
#' @export
mgeom <- function(prob) {
  c(mean     = (1 - prob) / prob,
    variance = (1 - prob) / prob^2)
}

#' @rdname disc.moments
#' @export
mnbinom <- function(size, prob) {
  c(mean     = size * (1 - prob) / prob,
    variance = size * (1 - prob) / prob^2)
}

#' @rdname disc.moments
#' @export
mhyper <- function(m, n, k) {
  N <- m + n
  c(mean     = k * m / N,
    variance = k * m / N * n / N * (N - k) / (N - 1))
}

#' @rdname disc.moments
#' @export
mbenford <- function(ndigits = 1) {
  d <- if (ndigits == 1) 1:9 else 10:99
  p <- log10(1 + 1/d)
  mu <- sum(d * p)
  c(mean     = mu,
    variance = sum(d^2 * p) - mu^2)
}