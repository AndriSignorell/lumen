
#' "Reverse" Gumbel Distribution Functions
#' 
#' Density, distribution function, quantile function and random generation for
#' the \dQuote{Reverse} Gumbel distribution with parameters \code{location} and
#' \code{scale}.
#' 
#' 
#' @name drevgumbel
#' @aliases dRevGumbel pRevGumbel qRevGumbel qRevGumbelExp rRevGumbel
#' @param x,q numeric vector of abscissa (or quantile) values at which to
#' evaluate the density or distribution function.
#' @param p numeric vector of probabilities at which to evaluate the quantile
#' function.
#' @param location location of the distribution
#' @param scale scale (\eqn{> 0}) of the distribution.
#' @param n number of random variates, i.e., \code{\link{length}} of resulting
#' vector of \code{rRevGumbel(..)}.
#' @return a numeric vector, of the same length as \code{x}, \code{q}, or
#' \code{p} for the first three functions, and of length \code{n} for
#' \code{rRevGumbel()}.
#' @seealso the \code{\link{Weibull}} distribution functions in \R's
#' \pkg{stats} package.
#' @note
#' Based on code by Werner Stahel; partly inspired by package \pkg{VGAM}. Martin
#' Maechler for numeric cosmetic.
#' 
#' @examples
#' 
#' curve(pRevGumbel(x, scale= 1/2), -3,2, n=1001, col=1, lwd=2,
#'       main = "RevGumbel(x, scale = 1/2)")
#' abline(h=0:1, v = 0, lty=3, col = "gray30")
#' curve(dRevGumbel(x, scale= 1/2),       n=1001, add=TRUE,
#'       col = (col.d <- adjustcolor(2, 0.5)), lwd=3)
#' legend("left", c("cdf","pdf"), col=c("black", col.d), lwd=2:3, bty="n")
#' 
#' med <- qRevGumbel(0.5, scale=1/2)
#' cat("The median is:",  format(med),"\n")
#' 

#' @rdname drevgumbel
#' @family dist.extreme
#' @concept distributions
#' @concept extreme-value-theory
#'
#'
#' @export
dRevGumbel <- function (x, location = 0, scale = 1) {
  # from VGAM  -- if (is.null(x)) FALSE else ifelse(is.na(x), FALSE, x)
  if (!isNumeric(scale, positive=TRUE))
    stop("\"scale\" must be positive")
  temp = exp((x - location)/scale)
  temp * exp(-temp)/scale
}

#' @rdname drevgumbel
#' @export
pRevGumbel <- function (q, location = 0, scale = 1) {
  
  if (!isNumeric(scale, positive=TRUE))
    stop("\"scale\" must be positive")
  1-exp(-exp((q - location)/scale))
}

#' @rdname drevgumbel
#' @export
qRevGumbel <- function (p, location = 0, scale = 1)
{
  if (!isNumeric(scale, positive=TRUE))
    stop("\"scale\" must be positive")
  location + scale * log(-log(p))
}


#' @rdname drevgumbel
#' @export
rRevGumbel <- function (n, location = 0, scale = 1)
{
  if (!isNumeric(scale, positive=TRUE, integer.valued=TRUE))
    stop("bad input for argument \"n\"")
  if (!isNumeric(scale, positive=TRUE))
    stop("\"scale\" must be positive")
  location + scale * log(-log(runif(n)))
}

#' @rdname drevgumbel
#' @export
qRevGumbelExp <- function (p) exp(qRevGumbel(p))

