
#' Confidence Interval for a Pearson Correlation
#'
#' Computes confidence intervals for a Pearson correlation coefficient
#' using Fisher's \eqn{z}-transformation. Since the sampling distribution
#' of the correlation coefficient is not normally distributed, Fisher's
#' transformation is applied to obtain approximately normally distributed
#' values, from which confidence intervals can be derived.
#'
#' @param rho numeric; Pearson correlation coefficient. Must be a single
#'   value in the interval \eqn{[-1, 1]}.
#' @param n integer; sample size used to estimate the correlation.
#'   Must be at least 3.
#' @param conf.level confidence level of the interval, a single number in
#'   \eqn{(0, 1)}. Defaults to \code{0.95}; unlike the estimators in this
#'   package it cannot be \code{NA}, since the interval is all this function
#'   computes.
#' @param sides character string specifying the sidedness of the confidence
#'   interval (one of \code{"two.sided"} (default), \code{"left"} or
#'   \code{"right"}). See \code{DescToolsX::ConfidenceIntervals}.
#'
#' @details
#' Fisher's \eqn{z}-transformation is defined as
#' \deqn{z = \frac{1}{2} \log\left(\frac{1 + r}{1 - r}\right),}
#' which stabilizes the variance of the correlation coefficient. The
#' transformed values are approximately normally distributed with standard
#' error \eqn{1 / \sqrt{n - 3}}. Confidence intervals are computed on the
#' transformed scale and then back-transformed.
#'
#' A correlation lies in \eqn{[-1, 1]}, so the open side of a one-sided
#' interval is reported at that boundary rather than at \eqn{\pm\infty}.
#'
#' The transformation has nothing to say in two situations, and both are
#' answered with \code{NA} bounds and a warning rather than with an interval
#' that only looks informative: at \eqn{|\rho| = 1}, where \eqn{z} is
#' infinite and the interval would collapse onto the estimate and thereby
#' rule out every other value, and at \eqn{n = 3}, where the standard error
#' \eqn{1/\sqrt{n-3}} is infinite and the interval would be the whole range.
#'
#' @section Argument name:
#' This function used to take \code{alternative} with the values
#' \code{"less"} and \code{"greater"}. It now takes \code{sides}, like every
#' other interval function in the package, and \code{sides} names the side
#' carrying the \emph{finite} bound - so \code{"left"} is the former
#' \code{"greater"} and \code{"right"} the former \code{"less"}. The values
#' produced are unchanged.
#'
#' @return A named numeric vector with elements:
#' \describe{
#'   \item{\code{est}}{point estimate (correlation coefficient \code{rho}).}
#'   \item{\code{lci}}{lower confidence interval bound.}
#'   \item{\code{uci}}{upper confidence interval bound.}
#' }
#'
#' @note Based on code by William Revelle, adapted to conform to package standards.
#'
#' @examples
#' # Confidence interval for a single correlation
#' corCI(0.5, n = 30)
#'
#' # one-sided: "left" carries the finite lower bound
#' corCI(0.5, n = 30, sides = "left")
#'
#' # Compare multiple correlations
#' r <- seq(0, 0.9, by = 0.1)
#' t(sapply(r, corCI, n = 30))
#'
#' @seealso \code{\link{fisherZ}}, \code{\link{fisherZInv}}, \code{\link{cor.test}}
#'
#' @family ci.correlation  
#' @concept correlation  
#' @concept confidence-interval
#'
#'
#' @export
corCI <- function(rho, n, conf.level = 0.95,
                  sides = c("two.sided", "left", "right")) {
  
  sides <- match.arg(sides)
  
  if (!is.numeric(rho) || length(rho) != 1L || !is.finite(rho))
    stop("'rho' must be a single finite numeric value")
  
  if (abs(rho) > 1)
    stop("'rho' must be between -1 and 1")
  
  # 'n' was compared against 3 without ever being checked, so corCI(0.5, "a")
  # failed on the comparison rather than on the argument
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) || n %% 1 != 0)
    stop("'n' must be a single whole number")
  
  if (n < 3L)
    stop("'n' must be at least 3")
  
  conf.level <- checkConfLevel(conf.level)
  
  # The shared check admits NA, because most functions in the suite use it
  # to mean "no interval wanted". Here there is nothing else to return.
  if (is.na(conf.level))
    stop("'conf.level' must be a number in (0, 1); corCI() computes nothing ",
         "but the interval")
  
  if (sides != "two.sided" && conf.level <= 0.5)
    stop(gettextf(
      "a one-sided interval needs 'conf.level' above 0.5, not %g",
      conf.level), domain = NA)
  
  # Two cases where the z-transformation is undefined. Both used to be
  # answered with something that looks like an interval: (rho, rho) at
  # |rho| = 1, which rules out every value below 1 and is a claim no sample
  # supports, and the whole range at n = 3, where 1/sqrt(n-3) is infinite.
  # cramerV(), spearmanCor() and kappaM() report NA in the same situations.
  if (abs(rho) >= 1) {
    warning("the z-transformation cannot bound a perfect correlation; ",
            "no interval computed", call. = FALSE)
    return(c(est = rho, lci = NA_real_, uci = NA_real_))
  }
  
  if (n == 3L) {
    warning("the z-transformation needs more than 3 observations; ",
            "no interval computed", call. = FALSE)
    return(c(est = rho, lci = NA_real_, uci = NA_real_))
  }
  
  # A one-sided bound at level gamma is the corresponding end of the
  # two-sided interval at level 2*gamma - 1, which is how every other
  # function in the suite writes it. The number is the qnorm(conf.level)
  # the two one-sided branches used to compute separately.
  confAdj <- if (sides == "two.sided") conf.level else 2 * conf.level - 1
  
  z     <- fisherZ(rho)
  sigma <- 1 / sqrt(n - 3)
  
  ci <- fisherZInv(z + c(-1, 1) * sigma * qnorm(1 - (1 - confAdj) / 2))
  
  # a correlation lives in [-1, 1]; applySides() clamps to that range and
  # closes the open side there. The former branches put +/-Inf on the
  # z-scale and let fisherZInv() map it to +/-1 - same numbers, but only
  # by way of tanh() rather than by saying so.
  c(est = rho, applySides(ci, sides, lo = -1, hi = 1))
  
}
