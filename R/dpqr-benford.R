#' Benford's Distribution
#' 
#' Density, distribution function, quantile function, and random generation for
#' Benford's distribution.
#' 
#' Benford's distribution (the significant-digit law) describes the 
#' probability distribution of the leading significant digit(s) in 
#' many naturally occurring numerical datasets. The probability of 
#' the first digit being d follows log10(1 + 1/d), with support on 
#' \{1, ..., 9\} for one leading digit or \{10, ..., 99\} for two leading 
#' digits. It is commonly applied in fraud detection and data quality 
#' assessment.
#' 
#' Benford's Law (aka \emph{the significant-digit law}) is the empirical
#' observation that in many naturally occurring tables of numerical data, the
#' leading significant (nonzero) digit is not uniformly distributed in
#' \eqn{\{1,2,\ldots,9\}}{1:9}. Instead, the leading significant digit
#' (\eqn{=D}, say) obeys the law 
#' \deqn{P(D=d) = \log_{10} \left( 1 + \frac1d \right)}{ P(D=d) = log10(1 + 1/d)} 
#' for \eqn{d=1,\ldots,9}. This means the
#' probability the first significant digit is 1 is approximately \eqn{0.301},
#' etc.
#' 
#' Benford's Law was apparently first discovered in 1881 by
#' astronomer/mathematician S. Newcomb. It started with the observation that the
#' pages of a book of logarithms were dirtiest at the beginning and
#' progressively cleaner throughout. In 1938, a General Electric physicist
#' called F. Benford rediscovered the law on this same observation. Over
#' several years he collected data from sources as different as
#' atomic weights, baseball statistics, numerical data from \emph{Reader's
#' Digest}, and drainage areas of rivers.
#' 
#' Applications of Benford's Law have been as diverse as fraud detection
#' in accounting and the design of computers.
#' 
#' @name dpqr-benford
#' @aliases Benford benford dbenford pbenford qbenford rbenford
#' 
#' @param x,q a vector of quantiles, see \code{nDigits}
#' 
#' @param p a vector of probabilities
#' @param n number of observations. A single positive integer.  Else if
#' \code{length(n) > 1} then the length is taken to be the number required.
#' 
#' @param nDigits number of leading digits, either 1 or 2.  If 1 then the
#' support of the distribution is \{1, ..., 9\}, else \{10, ..., 99\}.
#' 
#' @param log logical; if \code{TRUE}, densities are given as \code{log(d)}
#' 
#' @param log.p logical; if \code{TRUE}, probabilities \code{p} are given as
#' \code{log(p)}
#' 
#' @param lower.tail logical; if \code{TRUE} (default), probabilities are
#' \eqn{P(X \le x)}, otherwise \eqn{P(X > x)}
#' 
#' @return \code{dbenford()} gives the density, \code{pbenford()} gives the
#' distribution function, \code{qbenford()} gives the quantile function, and
#' \code{rbenford()} generates random deviates.
#' 
#' @references
#' Benford, F. (1938) The Law of Anomalous Numbers. \emph{Proceedings of the
#' American Philosophical Society}, \bold{78}, 551--572.
#' 
#' Newcomb, S. (1881) Note on the Frequency of Use of the Different Digits in
#' Natural Numbers. \emph{American Journal of Mathematics}, \bold{4}, 39--40.
#' 
#' @note
#' Based on code by T. W. Yee previously published in
#' the \pkg{VGAM} package, adapted to conform to package standards. 
#' 
#' @examples
#' dbenford(x <- c(0:10, NA, NaN, -Inf, Inf))
#' pbenford(x)
#' 
#' xx <- 1:9
#' barplot(dbenford(xx), col = "lightblue", las = 1, xlab = "Leading digit",
#'         ylab = "Probability", names.arg = as.character(xx),
#'         main = "Benford's distribution")
#' 
#' hist(rbenford(n = 1000), border = "blue", prob = TRUE,
#'      main = "1000 random variates from Benford's distribution",
#'      xlab = "Leading digit", sub = "Red is the true probability",
#'      breaks = 0:9 + 0.5, ylim = c(0, 0.35), xlim = c(0, 10.0))
#' lines(xx, dbenford(xx), col = "red", type = "h")
#' points(xx, dbenford(xx), col = "red")
#' 
#' @seealso [distributions-overview]
#' @concept distribution-function
#' @concept goodness-of-fit
#'
#' @rdname dpqr-benford
#' @export
dbenford <- function(x, nDigits = 1, log = FALSE) {
  if (!isNumeric(nDigits, isPositive = TRUE, isIntegerValued = TRUE) ||
      nDigits > 2)
    stop("argument 'nDigits' must be 1 or 2")
  lowerlimit <- ifelse(nDigits == 1, 1, 10)
  upperlimit <- ifelse(nDigits == 1, 9, 99)

  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)

  ans <- x * NA
  indexTF <- is.finite(x) & (x >= lowerlimit)

  ans[indexTF] <- log10(1 + 1/x[indexTF])
  ans[!is.na(x) & !is.nan(x) &
        ((x < lowerlimit) |
           (x > upperlimit) |
           (x != round(x)))] <- 0.0
  if (log.arg) log(ans) else ans
}


#' @rdname dpqr-benford
#' @export
pbenford <- function(q, nDigits = 1, lower.tail = TRUE, log.p = FALSE) {
  if (!isNumeric(nDigits, isPositive = TRUE, isIntegerValued = TRUE) ||
      nDigits > 2)
    stop("argument 'nDigits' must be 1 or 2")
  if (!is.logical(lower.tail) || length(lower.tail) != 1)
    stop("bad input for argument 'lower.tail'")
  if (!is.logical(log.p) || length(log.p) != 1)
    stop("bad input for argument 'log.p'")

  lowerlimit <- ifelse(nDigits == 1, 1, 10)
  upperlimit <- ifelse(nDigits == 1, 9, 99)

  ans <- q * NA
  floorq <- floor(q)
  indexTF <- is.finite(q) & (floorq >= lowerlimit)
  ans[indexTF] <- log10(1 + floorq[indexTF]) -
    ifelse(nDigits == 1, 0, 1)
  ans[!is.na(q) & !is.nan(q) & (q >= upperlimit)] <- 1
  ans[!is.na(q) & !is.nan(q) & (q <  lowerlimit)] <- 0

  if (!lower.tail) ans <- 1 - ans
  if (log.p) log(ans) else ans
}


#' @rdname dpqr-benford
#' @export
qbenford <- function(p, nDigits = 1, lower.tail = TRUE, log.p = FALSE) {
  if (!isNumeric(nDigits, isPositive = TRUE, isIntegerValued = TRUE) ||
      nDigits > 2)
    stop("argument 'nDigits' must be 1 or 2")
  if (!is.logical(lower.tail) || length(lower.tail) != 1)
    stop("bad input for argument 'lower.tail'")
  if (!is.logical(log.p) || length(log.p) != 1)
    stop("bad input for argument 'log.p'")

  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p

  lowerlimit <- ifelse(nDigits == 1, 1, 10)
  upperlimit <- ifelse(nDigits == 1, 9, 99)
  bad <- !is.na(p) & !is.nan(p) & ((p < 0) | (p > 1))
  if (any(bad))
    stop("bad input for argument 'p'")

  ans <- rep(lowerlimit, length = length(p))
  for (ii in (lowerlimit+1):upperlimit) {
    indexTF <- is.finite(p) &
      (pbenford(ii-1, nDigits = nDigits) < p) &
      (p <= pbenford(ii, nDigits = nDigits))
    ans[indexTF] <- ii
  }

  ans[ is.na(p) |  is.nan(p)] <- NA
  ans[!is.na(p) & !is.nan(p) & (p == 0)] <- lowerlimit
  ans[!is.na(p) & !is.nan(p) & (p == 1)] <- upperlimit
  ans
}


#' @rdname dpqr-benford
#' @export
rbenford <- function(n, nDigits = 1) {
  if (!isNumeric(nDigits, isPositive = TRUE, isIntegerValued = TRUE) ||
      nDigits > 2)
    stop("argument 'nDigits' must be 1 or 2")
  lowerlimit <- ifelse(nDigits == 1, 1, 10)
  upperlimit <- ifelse(nDigits == 1, 9, 99)
  use.n <- if ((length.n <- length(n)) > 1) length.n else
    if (!isNumeric(n, isIntegerValued = TRUE,
                   isPositive = TRUE))
      stop("bad input for argument 'n'") else n
  myrunif <- runif(use.n)

  ans <- rep(lowerlimit, length = use.n)
  for (ii in (lowerlimit+1):upperlimit) {
    indexTF <- (pbenford(ii-1, nDigits = nDigits) < myrunif) &
      (myrunif <= pbenford(ii, nDigits = nDigits))
    ans[indexTF] <- ii
  }
  ans
}
