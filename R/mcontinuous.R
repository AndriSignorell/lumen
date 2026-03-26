
#' Mean and Variance of Continuous Distributions
#'
#' Computes the mean and variance of common continuous distributions given
#' their parameters.
#'
#' @param mean mean of the normal distribution.
#' @param sd standard deviation of the normal distribution.
#' @param rate rate parameter (1/mean) of the exponential distribution.
#' @param shape, rate shape and scale parameters of the gamma and beta 
#'   distributions.
#' @param shape1,shape2 shape parameters of the beta distribution 
#'   (\eqn{\alpha} and \eqn{\beta}).
#' @param meanlog,sdlog mean and standard deviation on the log scale 
#'   (log-normal distribution).
#' @param df degrees of freedom (chi-squared and t-distribution).
#' @param df1,df2 numerator and denominator degrees of freedom 
#'   (F-distribution).
#' @param min lower limit of the triangular distribution.
#' @param max upper limit of the triangular distribution.
#' @param mode mode of the triangular distribution.
#'
#' @return A named numeric vector with elements \code{mean} and
#'   \code{variance}. Returns \code{NA} where moments do not exist.
#'
#' @details
#' \strong{Normal:}\cr
#' \eqn{\mu = \mu}\cr
#' \eqn{\mathrm{Var}(X) = \sigma^2}
#'
#' \strong{Exponential:}\cr
#' \eqn{\mu = \frac{1}{\lambda}}\cr
#' \eqn{\mathrm{Var}(X) = \frac{1}{\lambda^2}}
#'
#' \strong{Gamma:}\cr
#' \eqn{\mu = \frac{\alpha}{\beta}}\cr
#' \eqn{\mathrm{Var}(X) = \frac{\alpha}{\beta^2}}
#'
#' \strong{Log-Normal:}\cr
#' \eqn{\mu = \exp(\mu_{log} + \frac{1}{2}\sigma_{log}^2)}\cr
#' \eqn{\mathrm{Var}(X) = (\exp(\sigma_{log}^2) - 1) \cdot \exp(2\mu_{log} + \sigma_{log}^2)}
#'
#' \strong{Beta:}\cr
#' \eqn{\mu = \frac{\alpha}{\alpha + \beta}}\cr
#' \eqn{\mathrm{Var}(X) = \frac{\alpha\beta}{(\alpha+\beta)^2(\alpha+\beta+1)}}
#'
#' \strong{Chi-Squared:}\cr
#' \eqn{\mu = \nu}\cr
#' \eqn{\mathrm{Var}(X) = 2\nu}
#'
#' \strong{t-Distribution:}\cr
#' \eqn{\mu = 0 \quad (\nu > 1)}\cr
#' \eqn{\mathrm{Var}(X) = \frac{\nu}{\nu-2} \quad (\nu > 2)}
#'
#' \strong{F-Distribution:}\cr
#' \eqn{\mu = \frac{n_2}{n_2-2} \quad (n_2 > 2)}\cr
#' \eqn{\mathrm{Var}(X) = \frac{2n_2^2(n_1+n_2-2)}{n_1(n_2-2)^2(n_2-4)} \quad (n_2 > 4)}
#'
#' \strong{Triangular:}\cr
#' \eqn{\mu = \frac{a + b + c}{3}}\cr
#' \eqn{\mathrm{Var}(X) = \frac{a^2 + b^2 + c^2 - ab - ac - bc}{18}}\cr
#' where \eqn{a} = \code{min}, \eqn{b} = \code{max} and \eqn{c} = \code{mode}.
#'
#' @seealso \code{\link[stats]{dnorm}}, \code{\link[stats]{dexp}},
#'   \code{\link[stats]{dgamma}}, \code{\link[stats]{dlnorm}},
#'   \code{\link[stats]{dbeta}}, \code{\link[stats]{dchisq}},
#'   \code{\link[stats]{dt}}, \code{\link[stats]{df}}
#'
#' @references
#' Casella, G. and Berger, R. L. (2002) \emph{Statistical Inference}.
#' Duxbury.
#'
#' Johnson, N. L., Kotz, S. and Balakrishnan, N. (1994)
#' \emph{Continuous Univariate Distributions}, Vol. 1. Wiley.
#'
#' Johnson, N. L., Kotz, S. and Balakrishnan, N. (1995)
#' \emph{Continuous Univariate Distributions}, Vol. 2. Wiley.
#'
#' @family topic.distributions.moments
#' @concept continuous distribution
#' @concept moments
#'
#' @examples
#' mnorm(mean = 0, sd = 1)
#' mexp(rate = 0.5)
#' mgamma(shape = 2, rate = 0.5)
#' mlnorm(meanlog = 0, sdlog = 1)
#' mbeta(shape1 = 2, shape2 = 3)
#' mchisq(df = 4)
#' mt(df = 5)
#' mf(df1 = 5, df2 = 10)
#' mtri(min = 0, max = 1, mode = 0.5)
#' mtri(min = 2, max = 10, mode = 4)
#'
#' @name cont.moments
NULL

#' @rdname cont.moments
#' @export
mnorm <- function(mean, sd) {
  c(mean     = mean,
    variance = sd^2)
}

#' @rdname cont.moments
#' @export
mexp <- function(rate) {
  c(mean     = 1 / rate,
    variance = 1 / rate^2)
}

#' @rdname cont.moments
#' @export
mgamma <- function(shape, rate) {
  c(mean     = shape / rate,
    variance = shape / rate^2)
}

#' @rdname cont.moments
#' @export
mlnorm <- function(meanlog, sdlog) {
  c(mean     = exp(meanlog + 0.5 * sdlog^2),
    variance = (exp(sdlog^2) - 1) * exp(2 * meanlog + sdlog^2))
}

#' @rdname cont.moments
#' @export
mbeta <- function(shape1, shape2) {
  c(mean     = shape1 / (shape1 + shape2),
    variance = (shape1 * shape2) /
      ((shape1 + shape2)^2 * (shape1 + shape2 + 1)))
}

#' @rdname cont.moments
#' @export
mchisq <- function(df) {
  c(mean     = df,
    variance = 2 * df)
}

#' @rdname cont.moments
#' @export
mt <- function(df) {
  c(mean     = if (df > 1) 0 else NA_real_,
    variance = if (df > 2) df / (df - 2) else NA_real_)
}

#' @rdname cont.moments
#' @export
mf <- function(df1, df2) {
  c(mean     = if (df2 > 2) df2 / (df2 - 2) else NA_real_,
    variance = if (df2 > 4)
      2 * df2^2 * (df1 + df2 - 2) /
      (df1 * (df2 - 2)^2 * (df2 - 4))
    else NA_real_)
}

#' @rdname cont.moments
#' @export
mtri <- function(min = 0, max = 1, mode = 0.5) {
  c(mean     = (min + max + mode) / 3,
    variance = (min^2 + max^2 + mode^2 - 
                  min*max - min*mode - max*mode) / 18)
}

