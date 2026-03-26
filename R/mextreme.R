
#' Mean and Variance of Extreme Value Distributions
#'
#' Computes the mean and variance of common extreme value distributions
#' given their parameters.
#'
#' @param loc location parameter.
#' @param scale scale parameter.
#' @param shape shape parameter.
#'
#' @return A named numeric vector with elements \code{mean} and
#'   \code{variance}. Returns \code{NA} where moments do not exist.
#'
#' @details
#' \strong{Gumbel:}\cr
#' \eqn{\mu = a + b\gamma}\cr
#' \eqn{\mathrm{Var}(X) = \frac{\pi^2}{6} b^2}\cr
#' where \eqn{\gamma \approx 0.5772} is the Euler-Mascheroni constant,
#' \eqn{a} = \code{loc} and \eqn{b} = \code{scale}.
#'
#' \strong{Fréchet:}\cr
#' \eqn{\mu = a + b\,\Gamma(1 - 1/s) \quad (s > 1)}\cr
#' \eqn{\mathrm{Var}(X) = b^2\left[\Gamma(1 - 2/s) - 
#'   \Gamma(1 - 1/s)^2\right] \quad (s > 2)}\cr
#' where \eqn{a} = \code{loc}, \eqn{b} = \code{scale} and 
#' \eqn{s} = \code{shape}.
#'
#' \strong{Reverse Weibull:}\cr
#' \eqn{\mu = a - b\,\Gamma(1 + 1/s)}\cr
#' \eqn{\mathrm{Var}(X) = b^2\left[\Gamma(1 + 2/s) - 
#'   \Gamma(1 + 1/s)^2\right]}\cr
#' where \eqn{a} = \code{loc}, \eqn{b} = \code{scale} and 
#' \eqn{s} = \code{shape}.
#'
#' \strong{GEV:}\cr
#' For \eqn{s = 0} (Gumbel), \eqn{s > 0} (Fréchet) and \eqn{s < 0}
#' (Reverse Weibull), the moments are computed accordingly. Mean exists
#' for \eqn{s < 1}, variance for \eqn{s < 1/2}.
#'
#' \strong{GPD:}\cr
#' \eqn{\mu = a + \frac{b}{1-s} \quad (s < 1)}\cr
#' \eqn{\mathrm{Var}(X) = \frac{b^2}{(1-s)^2(1-2s)} \quad (s < 1/2)}\cr
#' where \eqn{a} = \code{loc}, \eqn{b} = \code{scale} and 
#' \eqn{s} = \code{shape}.
#'
#' \strong{Gompertz:}\cr
#' For \eqn{a > 0} the moments are computed numerically via integration.
#' For \eqn{a = 0} the distribution reduces to the exponential with
#' \eqn{\mu = 1/b} and \eqn{\mathrm{Var}(X) = 1/b^2}.
#' For \eqn{a < 0} the moments do not exist (\code{NA}).
#' 
#' @seealso \code{\link{dgumbel}}, \code{\link{dfrechet}},
#'   \code{\link{drweibull}}, \code{\link{dgev}}, \code{\link{dgpd}}
#'
#' @references
#' Coles, S. (2001) \emph{An Introduction to Statistical Modeling of
#' Extreme Values}. Springer.
#'
#' Kotz, S. and Nadarajah, S. (2000) \emph{Extreme Value Distributions}.
#' Imperial College Press.
#'
#' @family topic.distributions.moments
#' @concept continuous distribution
#' @concept extreme value theory
#' @concept moments
#'
#' @examples
#' mgumbel(loc = 0, scale = 1)
#' mfrechet(loc = 0, scale = 1, shape = 3)
#' mrweibull(loc = 0, scale = 1, shape = 2)
#' mgev(loc = 0, scale = 1, shape = 0)
#' mgev(loc = 0, scale = 1, shape = 0.3)
#' mgev(loc = 0, scale = 1, shape = -0.3)
#' mgpd(loc = 0, scale = 1, shape = 0.3)
#'
#' @name ev.moments
NULL

#' @rdname ev.moments
#' @export
mgumbel <- function(loc = 0, scale = 1) {
  c(mean     = loc + scale * 0.5772156649015329,
    variance = pi^2 / 6 * scale^2)
}

#' @rdname ev.moments
#' @export
mfrechet <- function(loc = 0, scale = 1, shape = 1) {
  c(mean     = if (shape > 1)
    loc + scale * gamma(1 - 1/shape)
    else NA_real_,
    variance = if (shape > 2)
      scale^2 * (gamma(1 - 2/shape) - gamma(1 - 1/shape)^2)
    else NA_real_)
}

#' @rdname ev.moments
#' @export
mrweibull <- function(loc = 0, scale = 1, shape = 1) {
  c(mean     = loc - scale * gamma(1 + 1/shape),
    variance = scale^2 * (gamma(1 + 2/shape) - gamma(1 + 1/shape)^2))
}

#' @rdname ev.moments
#' @export
mgev <- function(loc = 0, scale = 1, shape = 0) {
  if (shape == 0) {
    # Gumbel
    c(mean     = loc + scale * 0.5772156649015329,
      variance = pi^2 / 6 * scale^2)
  } else if (shape < 1 && shape != 0) {
    g1 <- gamma(1 - shape)
    g2 <- gamma(1 - 2 * shape)
    c(mean     = loc + scale * (g1 - 1) / shape,
      variance = if (shape < 0.5)
        scale^2 * (g2 - g1^2) / shape^2
      else NA_real_)
  } else {
    c(mean = NA_real_, variance = NA_real_)
  }
}

#' @rdname ev.moments
#' @export
mgpd <- function(loc = 0, scale = 1, shape = 0) {
  c(mean     = if (shape < 1)
    loc + scale / (1 - shape)
    else NA_real_,
    variance = if (shape < 0.5)
      scale^2 / ((1 - shape)^2 * (1 - 2 * shape))
    else NA_real_)
}



#' @rdname ev.moments
#' @export
mgompertz <- function(shape, rate = 1) {
  if (shape == 0) {
    # exponential special case
    c(mean     = 1 / rate,
      variance = 1 / rate^2)
  } else if (shape > 0) {
    mu <- integrate(function(x)
      pgompertz(x, shape = shape, rate = rate,
                lower.tail = FALSE),
      0, Inf)$value
    c(mean     = mu,
      variance = integrate(function(x)
        2 * x * pgompertz(x, shape = shape, rate = rate,
                          lower.tail = FALSE),
        0, Inf)$value - mu^2)
  } else {
    # shape < 0: non-zero probability of immortality, moments do not exist
    c(mean = NA_real_, variance = NA_real_)
  }
}

