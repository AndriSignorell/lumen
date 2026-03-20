
#' The Gumbel Distribution
#' 
#' Density function, distribution function, quantile function and random
#' generation for the Gumbel distribution with location and scale parameters.
#' 
#' The Gumbel distribution function with parameters \eqn{\code{loc} = a} and
#' \eqn{\code{scale} = b} is 
#' \deqn{G(z) = \exp\left\{-\exp\left[-\left(\frac{z-a}{b}\right)\right]\right\}}{G(z) = exp(-exp(-(z-a)/b))}
#' for all real \eqn{z}, where \eqn{b > 0}.
#'  
#' @name dgumbel
#' @aliases dgumbel pgumbel qgumbel rgumbel
#' @param x,q Vector of quantiles.
#' @param p Vector of probabilities.
#' @param n Number of observations.
#' @param loc,scale Location and scale parameters (can be given as vectors).
#' @param log Logical; if \code{TRUE}, the log density is returned.
#' @param lower.tail Logical; if \code{TRUE} (default), probabilities are 
#' \verb{P[X <= x]}, otherwise, P\verb{[X > x]}
#' @return \code{dgumbel} gives the density function, \code{pgumbel} gives the
#' distribution function, \code{qgumbel} gives the quantile function, and
#' \code{rgumbel} generates random deviates.
#' 
#' @family topic.distributions
#' @concept continuous distribution
#' @concept extreme value theory
#' @concept Gumbel
#' @concept dpqr
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


#' @rdname dgumbel
#' @export
"dgumbel"<-
  function(x, loc = 0, scale = 1, log = FALSE)
  {
    dgev(x, loc = loc, scale = scale, shape = 0, log = log)
  }


#' @rdname dgumbel
#' @export
"pgumbel"<-
  function(q, loc = 0, scale = 1, lower.tail = TRUE)
  {
    pgev(q, loc = loc, scale = scale, shape = 0, lower.tail = lower.tail)
  }


#' @rdname dgumbel
#' @export
"qgumbel"<-
  function(p, loc = 0, scale = 1, lower.tail = TRUE)
  {
    qgev(p, loc = loc, scale = scale, shape = 0, lower.tail = lower.tail)
  }


#' @rdname dgumbel
#' @export
"rgumbel"<-
  function(n, loc = 0, scale = 1)
  {
    rgev(n, loc = loc, scale = scale, shape = 0)
  }




