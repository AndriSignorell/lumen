
#' Mean and Variance of Extreme Value Distributions
#'
#' Computes the mean and variance of common extreme value distributions
#' given their parameters.
#'
#' @param loc location parameter.
#' @param scale scale parameter.
#' @param shape,rate shape and rate parameters.
#'
#' @return A named numeric vector with elements `mean` and
#'   `variance`. Returns `NA` where moments do not exist.
#'
#' @details
#' \tabular{lll}{
#'   **Distribution** \tab **Mean** 
#'   \tab **Variance** \cr
#'   Gumbel \tab
#'     \eqn{a + b\gamma} \tab
#'     \eqn{\frac{\pi^2}{6}b^2} \cr
#'   Reverse Gumbel \tab
#'     \eqn{a - b\gamma} \tab
#'     \eqn{\frac{\pi^2}{6}b^2} \cr
#'   Fréchet \tab
#'     \eqn{a + b\Gamma(1 - 1/s) \quad (s > 1)} \tab
#'     \eqn{b^2\left[\Gamma(1 - 2/s) -
#'       \Gamma(1 - 1/s)^2\right] \quad (s > 2)} \cr
#'   Reverse Weibull \verb{  } \tab
#'     \eqn{a - b\Gamma(1 + 1/s)} \tab
#'     \eqn{b^2\left[\Gamma(1 + 2/s) -
#'       \Gamma(1 + 1/s)^2\right]} \cr
#'   GEV \tab
#'     \eqn{a + b\frac{\Gamma(1-s)-1}{s}
#'       \quad (s \ne 0,\ s < 1)} \verb{  } \tab
#'     \eqn{b^2\frac{\Gamma(1-2s)-\Gamma(1-s)^2}{s^2}
#'       \quad (s \ne 0,\ s < 1/2)} \cr
#'   GPD \tab
#'     \eqn{a + \frac{b}{1-s} \quad (s < 1)} \tab
#'     \eqn{\frac{b^2}{(1-s)^2(1-2s)} \quad (s < 1/2)} \cr
#'   Gompertz \tab
#'     numerical integration for \eqn{\alpha > 0} \tab dito \cr
#'     \tab  \eqn{1/\beta}  for \eqn{\alpha = 0}; 
#'     \tab \eqn{1/\beta^2} for \eqn{\alpha = 0}; \cr
#'     \tab `NA` for \eqn{\alpha < 0} \tab dito \cr
#' }
#'
#' For the first six distributions, \eqn{a} = `loc`,
#' \eqn{b} = `scale`, and \eqn{s} = `shape`. For the GEV with
#' \eqn{s = 0}, the Gumbel moments apply. Furthermore,
#' \eqn{\gamma \approx 0.5772} is the Euler-Mascheroni constant. For the
#' Gompertz distribution, \eqn{\alpha} = `shape` and
#' \eqn{\beta} = `rate`; moments for \eqn{\alpha > 0} are computed
#' numerically by integration.
#' 
#' @seealso [dgumbel()], [drevgumbel()], [dfrechet()],
#'   [drevweibull()], [dgev()], [dgpd()],
#'   [distributions-overview]
#'
#' @references
#' Coles, S. (2001) *An Introduction to Statistical Modeling of
#' Extreme Values*. Springer.
#'
#' Kotz, S. and Nadarajah, S. (2000) *Extreme Value Distributions*.
#' Imperial College Press.
#'
#' @concept distribution-summary
#' @concept extreme-value
#' @concept moment
#'
#' @examples
#' mgumbel(loc = 0, scale = 1)
#' mrevgumbel(loc = 0, scale = 1)
#' mfrechet(loc = 0, scale = 1, shape = 3)
#' mrevweibull(loc = 0, scale = 1, shape = 2)
#' mgev(loc = 0, scale = 1, shape = 0)
#' mgev(loc = 0, scale = 1, shape = 0.3)
#' mgev(loc = 0, scale = 1, shape = -0.3)
#' mgpd(loc = 0, scale = 1, shape = 0.3)
#'
#' @name extreme-value-moments
NULL

#' @rdname extreme-value-moments
#' @export
mgumbel <- function(loc = 0, scale = 1) {

  .assertScalar(loc)
  .assertScalar(scale, lower = 0, strictLower = TRUE)

  c(mean     = loc + scale * 0.5772156649015329,
    variance = pi^2 / 6 * scale^2)
}

#' @rdname extreme-value-moments
#' @export
mrevgumbel <- function(loc = 0, scale = 1) {

  .assertScalar(loc)
  .assertScalar(scale, lower = 0, strictLower = TRUE)

  # the reflection X = a - bY of the standard Gumbel Y flips the sign of the
  # scale term in the mean and leaves the variance untouched
  c(mean     = loc - scale * 0.5772156649015329,
    variance = pi^2 / 6 * scale^2)
}

#' @rdname extreme-value-moments
#' @export
mfrechet <- function(loc = 0, scale = 1, shape = 1) {

  .assertScalar(loc)
  .assertScalar(scale, lower = 0, strictLower = TRUE)
  .assertScalar(shape, lower = 0, strictLower = TRUE)

  c(mean     = if (shape > 1)
    loc + scale * gamma(1 - 1/shape)
    else NA_real_,
    variance = if (shape > 2)
      scale^2 * (gamma(1 - 2/shape) - gamma(1 - 1/shape)^2)
    else NA_real_)
}

#' @rdname extreme-value-moments
#' @export
mrevweibull <- function(loc = 0, scale = 1, shape = 1) {

  .assertScalar(loc)
  .assertScalar(scale, lower = 0, strictLower = TRUE)
  .assertScalar(shape, lower = 0, strictLower = TRUE)

  c(mean     = loc - scale * gamma(1 + 1/shape),
    variance = scale^2 * (gamma(1 + 2/shape) - gamma(1 + 1/shape)^2))
}

#' @rdname extreme-value-moments
#' @export
mgev <- function(loc = 0, scale = 1, shape = 0) {

  .assertScalar(loc)
  .assertScalar(scale, lower = 0, strictLower = TRUE)
  .assertScalar(shape)

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

#' @rdname extreme-value-moments
#' @export
mgpd <- function(loc = 0, scale = 1, shape = 0) {

  .assertScalar(loc)
  .assertScalar(scale, lower = 0, strictLower = TRUE)
  .assertScalar(shape)

  c(mean     = if (shape < 1)
    loc + scale / (1 - shape)
    else NA_real_,
    variance = if (shape < 0.5)
      scale^2 / ((1 - shape)^2 * (1 - 2 * shape))
    else NA_real_)
}



#' @rdname extreme-value-moments
#' @export
mgompertz <- function(shape, rate = 1) {

  .assertScalar(shape)
  .assertScalar(rate, lower = 0, strictLower = TRUE)
  
  if (shape == 0) {
    return(c(mean = 1 / rate,
             variance = 1 / rate^2))
    
  } else if (shape > 0) {
    
    mu <- integrate(function(x)
      pgompertz(x, shape = shape, rate = rate,
                lower.tail = FALSE),
      0, Inf)$value
    
    return(c(mean = mu,
             variance = integrate(function(x)
               2 * x * pgompertz(x, shape = shape, rate = rate,
                                 lower.tail = FALSE),
               0, Inf)$value - mu^2))
    
  } else {
    return(c(mean = NA_real_, variance = NA_real_))
  }
}
