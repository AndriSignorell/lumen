
#' Benford's Distribution
#' 
#' Density, distribution function, quantile function, and random generation for
#' Benford's distribution.
#' 
#' 
#' Benford's Law (aka \emph{the significant-digit law}) is the empirical
#' observation that in many naturally occuring tables of numerical data, the
#' leading significant (nonzero) digit is not uniformly distributed in
#' \eqn{\{1,2,\ldots,9\}}{1:9}. Instead, the leading significant digit
#' (\eqn{=D}, say) obeys the law 
#' \deqn{P(D=d) = \log_{10} \left( 1 + \frac1d \right)}{ P(D=d) = log10(1 + 1/d)} 
#' for \eqn{d=1,\ldots,9}. This means the
#' probability the first significant digit is 1 is approximately \eqn{0.301},
#' etc.
#' 
#' Benford's Law was apparently first discovered in 1881 by
#' astronomer/mathematician S. Newcombe. It started by the observation that the
#' pages of a book of logarithms were dirtiest at the beginning and
#' progressively cleaner throughout. In 1938, a General Electric physicist
#' called F. Benford rediscovered the law on this same observation. Over
#' several years he collected data from different sources as different as
#' atomic weights, baseball statistics, numerical data from \emph{Reader's
#' Digest}, and drainage areas of rivers.
#' 
#' Applications of Benford's Law has been as diverse as to the area of fraud
#' detection in accounting and the design computers.
#' 
#' @name dbenford
#' @aliases Benford dbenford pbenford qbenford rbenford
#' @param x,q Vector of quantiles.  See \code{ndigits}.
#' 
#' @param p vector of probabilities.
#' @param n number of observations. A single positive integer.  Else if
#' \code{length(n) > 1} then the length is taken to be the number required.
#' 
#' @param ndigits Number of leading digits, either 1 or 2.  If 1 then the
#' support of the distribution is \{1, ...,9\}, else \{10, ..., 99\}.
#' 
#' @param log,log.p Logical.  If \code{log.p = TRUE} then all probabilities
#' \code{p} are given as \code{log(p)}.
#' 
#' @return \code{dBenf} gives the density, \code{pBenf} gives the distribution
#' function, and \code{qBenf} gives the quantile function, and \code{rBenf}
#' generates random deviates.
#' 
#' @author T. W. Yee
#' 
#' @references
#' Benford, F. (1938) The Law of Anomalous Numbers. \emph{Proceedings of the
#' American Philosophical Society}, \bold{78}, 551--572.
#' 
#' Newcomb, S. (1881) Note on the Frequency of Use of the Different Digits in
#' Natural Numbers. \emph{American Journal of Mathematics}, \bold{4}, 39--40.
#' @source These functions were previously published as \code{dbenf()} etc. in
#' the \pkg{VGAM} package and have been integrated here without logical
#' changes.
#' 
#' @family topic.distributions
#' @concept discrete distribution
#' @concept Benford law
#' @concept significant digits
#' @concept fraud detection
#' @concept dpqr
#' 
#' @examples
#' 
#' dbenford(x <- c(0:10, NA, NaN, -Inf, Inf))
#' pbenford(x)
#' 
#' \dontrun{
#' xx <- 1:9
#' barplot(dbenford(xx), col = "lightblue", las = 1, xlab = "Leading digit",
#'         ylab = "Probability", names.arg = as.character(xx),
#'         main = paste("Benford's distribution",  sep = ""))
#' 
#' hist(rbenford(n = 1000), border = "blue", prob = TRUE,
#'      main = "1000 random variates from Benford's distribution",
#'      xlab = "Leading digit", sub="Red is the true probability",
#'      breaks = 0:9 + 0.5, ylim = c(0, 0.35), xlim = c(0, 10.0))
#' lines(xx, dbenford(xx), col = "red", type = "h")
#' points(xx, dbenford(xx), col = "red")
#' }
#' 

#' @rdname dbenford
#' @export
dbenford <- function(x, ndigits = 1, log = FALSE) {
  if (!isNumeric(ndigits, length.arg = 1,
                 positive = TRUE, integer.valued = TRUE) ||
      ndigits > 2)
    stop("argument 'ndigits' must be 1 or 2")
  lowerlimit <- ifelse(ndigits == 1, 1, 10)
  upperlimit <- ifelse(ndigits == 1, 9, 99)
  
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


#' @rdname dbenford
#' @export
rbenford <- function(n, ndigits = 1) {
  if (!isNumeric(ndigits, length.arg = 1,
                 positive = TRUE, integer.valued = TRUE) ||
      ndigits > 2)
    stop("argument 'ndigits' must be 1 or 2")
  lowerlimit <- ifelse(ndigits == 1, 1, 10)
  upperlimit <- ifelse(ndigits == 1, 9, 99)
  use.n <- if ((length.n <- length(n)) > 1) length.n else
    if (!isNumeric(n, integer.valued = TRUE,
                   length.arg = 1, positive = TRUE))
      stop("bad input for argument 'n'") else n
  myrunif <- runif(use.n)
  
  ans <- rep(lowerlimit, length = use.n)
  for (ii in (lowerlimit+1):upperlimit) {
    indexTF <- (pbenford(ii-1, ndigits = ndigits) < myrunif) &
      (myrunif <= pbenford(ii, ndigits = ndigits))
    ans[indexTF] <- ii
  }
  ans
}


#' @rdname dbenford
#' @export
pbenford <- function(q, ndigits = 1, log.p = FALSE) {
  if (!isNumeric(ndigits, length.arg = 1,
                 positive = TRUE, integer.valued = TRUE) ||
      ndigits > 2)
    stop("argument 'ndigits' must be 1 or 2")
  lowerlimit <- ifelse(ndigits == 1, 1, 10)
  upperlimit <- ifelse(ndigits == 1, 9, 99)
  
  ans <- q * NA
  floorq <- floor(q)
  indexTF <- is.finite(q) & (floorq >= lowerlimit)
  ans[indexTF] <- log10(1 + floorq[indexTF]) -
    ifelse(ndigits == 1, 0, 1)
  ans[!is.na(q) & !is.nan(q) & (q >= upperlimit)] <- 1
  ans[!is.na(q) & !is.nan(q) & (q <  lowerlimit)] <- 0
  if (log.p) log(ans) else ans
}


#' @rdname dbenford
#' @export
qbenford <- function(p, ndigits = 1) {
  if (!isNumeric(ndigits, length.arg = 1,
                 positive = TRUE, integer.valued = TRUE) ||
      ndigits > 2)
    stop("argument 'ndigits' must be 1 or 2")
  lowerlimit <- ifelse(ndigits == 1, 1, 10)
  upperlimit <- ifelse(ndigits == 1, 9, 99)
  bad <- !is.na(p) & !is.nan(p) & ((p < 0) | (p > 1))
  if (any(bad))
    stop("bad input for argument 'p'")
  
  ans <- rep(lowerlimit, length = length(p))
  for (ii in (lowerlimit+1):upperlimit) {
    indexTF <- is.finite(p) &
      (pbenford(ii-1, ndigits = ndigits) < p) &
      (p <= pbenford(ii, ndigits = ndigits))
    ans[indexTF] <- ii
  }
  
  ans[ is.na(p) |  is.nan(p)] <- NA
  ans[!is.na(p) & !is.nan(p) & (p == 0)] <- lowerlimit
  ans[!is.na(p) & !is.nan(p) & (p == 1)] <- upperlimit
  ans
}


