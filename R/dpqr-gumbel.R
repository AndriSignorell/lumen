
#' Gumbel Distribution
#' 
#' The Gumbel distribution, also known as the Type I extreme value 
#' distribution, is a continuous distribution used to model the maximum 
#' (or minimum) of a number of samples of various distributions. It is 
#' parameterized by location and scale and arises naturally in extreme 
#' value theory.
#' 
#' Density function, distribution function, quantile function and random
#' generation for the Gumbel distribution with location and scale parameters.
#' 
#' The Gumbel distribution function with parameters \eqn{\code{loc} = a} and
#' \eqn{\code{scale} = b} is 
#' \deqn{G(z) = \exp\left\{-\exp\left[-\left(\frac{z-a}{b}\right)\right]\right\}}{G(z) = exp(-exp(-(z-a)/b))}
#' for all real \eqn{z}, where \eqn{b > 0}.
#'  
#' @name dpqr-gumbel
#' @aliases dgumbel pgumbel qgumbel rgumbel
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations.
#' @param loc,scale location and scale parameters (can be given as vectors).
#' @param log logical; if \code{TRUE}, the log density is returned.
#' @param lower.tail logical; if \code{TRUE} (default), probabilities are 
#' \verb{P[X <= x]}, otherwise, P\verb{[X > x]}.
#' @return \code{dgumbel()} gives the density function, \code{pgumbel()} gives
#' the distribution function, \code{qgumbel()} gives the quantile function,
#' and \code{rgumbel()} generates random deviates.
#' 
#' @note
#' Based on code by Alec Stephenson previously published in
#' the \pkg{evd} package, adapted to conform to package standards.
#' 
#' @seealso [distributions-overview]
#' @concept distribution-function
#' @concept extreme-value
#' 
#' @examples
#' 
#' dgumbel(-1:2, -1, 0.5)
#' pgumbel(-1:2, -1, 0.5)
#' qgumbel(seq(0.9, 0.6, -0.1), 2, 0.5)
#' rgumbel(6, -1, 0.5)
#' p <- (1:9)/10
#' pgumbel(qgumbel(p, -1, 2), -1, 2)
#' ## [1] 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9
#' 


#' @rdname dpqr-gumbel
#' @export
dgumbel <- function(x, loc = 0, scale = 1, log = FALSE)
  {
    dgev(x, loc = loc, scale = scale, shape = 0, log = log)
  }


#' @rdname dpqr-gumbel
#' @export
pgumbel <- function(q, loc = 0, scale = 1, lower.tail = TRUE)
  {
    pgev(q, loc = loc, scale = scale, shape = 0, lower.tail = lower.tail)
  }


#' @rdname dpqr-gumbel
#' @export
qgumbel <- function(p, loc = 0, scale = 1, lower.tail = TRUE)
  {
    qgev(p, loc = loc, scale = scale, shape = 0, lower.tail = lower.tail)
  }


#' @rdname dpqr-gumbel
#' @export
rgumbel <- function(n, loc = 0, scale = 1)
  {
    rgev(n, loc = loc, scale = scale, shape = 0)
  }




