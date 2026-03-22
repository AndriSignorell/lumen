
#' Moses Test of Extreme Reactions
#' 
#' A nonparametric test comparing the spread (range) of two independent 
#' groups, assessing whether the control group shows greater variability 
#' than the treatment group.
#' 
#' Perform Moses test of extreme reactions, which is a distribution-free
#' non-parametric test for the difference between two independent groups in the
#' extremity of scores (in both directions) that the groups contain.  Scores
#' from both groups are pooled and converted to ranks, and the test statistic
#' is the span of scores (the range plus 1) in one of the groups chosen
#' arbitrarily. An exact probability is computed for the span and then
#' recomputed after dropping a specified number of extreme scores from each end
#' of its range. The exact one-tailed probability is calculated.
#' 
#' For two independent samples from a continuous field, this tests whether
#' extreme values are equally likely in both populations or if they are more
#' likely to occur in the population from which the sample with the larger
#' range was drawn.
#' 
#' Note that the ranks are calculated in decreasing mode.
#' 
#' @name mosesTest
#' @aliases mosesTest mosesTest.default mosesTest.formula
#' @param x numeric vector of data values. \code{x} will be treated as control
#' group. Non-finite (e.g. infinite or missing) values will be omitted.
#' @param y numeric vector of data values. \code{y} will be treated as
#' experiment group. Non-finite (e.g. infinite or missing) values will be
#' omitted.
#' @param formula a formula of the form \code{lhs ~ rhs} where \code{lhs} gives
#' the data values and rhs the corresponding groups.
#' @param data an optional matrix or data frame (or similar: see
#' \code{\link{model.frame}}) containing the variables in the formula
#' \code{formula}.  By default the variables are taken from
#' \code{environment(formula)}.
#' @param subset an optional vector specifying a subset of observations to be
#' used.
#' @param na.action a function which indicates what should happen when the data
#' contain \code{NA}s. Defaults to \code{getOption("na.action")}.
#' @param extreme integer, defines the number of extreme values to be dropped
#' from the control group before calculating the span. Default (\code{NULL}) is
#' the integer part of \code{0.05 * length(x)} or \code{1}, whichever is
#' greater. If extreme is too large, it will be cut down to
#' \code{floor(length(x)-2)/2}.
#' @param \dots further arguments to be passed to or from methods.
#' @return A list with class \dQuote{htest} containing the following
#' components: \item{statistic}{the value of the Moses Test statistic.}
#' \item{p.value}{the p-value for the test.} \item{method}{the character string
#' \dQuote{Moses Test of Extreme Reactions}.} \item{data.name}{a character
#' string giving the name(s) of the data.}
#' 
#' @seealso \code{\link{wilcox.test}}, \code{\link{ks.test}}
#' 
#' @references Moses, L.E. (1952) A Two-Sample Test, \emph{Psychometrika}, 17,
#' 239-247.
#' 
#' @family topic.nonparametricTests
#' @concept scale test
#' 
#' @examples
#' 
#' x <- c(0.80, 0.83, 1.89, 1.04, 1.45, 1.38, 1.91, 1.64, 0.73, 1.46)
#' y <- c(1.15, 0.88, 0.90, 0.74, 1.21)
#' 
#' mosesTest(x, y)
#' 
#' 
#' set.seed(1479)
#' x <- sample(1:20, 10, replace=TRUE)
#' y <- sample(5:25, 6, replace=TRUE)
#' 
#' mosesTest(x, y)
#' 


#' @rdname mosesTest
#' @export
mosesTest <- function (x, ...)  UseMethod("mosesTest")

# Extremreaktionen nach Moses: Nullhypothese: Die Spannweite der Werte ist
# in beiden Gruppen gleich gross. Die Werte beider Gruppen werden in eine gemeinsame
# Reihenfolge gebracht. Anschliessend werden ihnen Rangwerte zugewiesen.
# Eine der Gruppen (die Gruppe des Wertes, der in dem Dialogfeld
#                   Gruppen definieren als erstes angegeben ist) wird als Kontrollgruppe verwendet.
# Fuer diese Gruppe wird die Spannweite der Raenge als Differenz zwischen
# dem groessten und kleinsten Rangwert berechnet. Anhand dieser Spannweite errechnet
# sich die einseitige Signifikanz. Zusaetzlich wird der Test ein zweites
# Mal durchgefuehrt, wobei die Ausreisser der Gesamtstichprobe ausgeschlossen
# werden (insgesamt 5% der Faelle). Dabei werden sowohl die hoechsten als auch
# die niedrigsten Raenge entfernt. Das Testergebnis teilt die Anzahl der Faelle beider
# Gruppen, die Spannweiten und die entsprechenden einseitigen Signifikanzen
# fuer beide Tests (mit und ohne Ausreisser) mit. Fuer ein Beispiel siehe oben
# Abschnitt Moses-Test, S. 760.


#' @rdname mosesTest
#' @export
mosesTest.formula <- function (formula, data, subset, na.action, ...) {
  
  # this is a taken analogue to wilcox.test.formula
  
  if (missing(formula) || (length(formula) != 3L) || (length(attr(terms(formula[-2L]),
                                                                  "term.labels")) != 1L))
    stop("'formula' missing or incorrect")
  m <- match.call(expand.dots = FALSE)
  if (is.matrix(eval(m$data, parent.frame())))
    m$data <- as.data.frame(data)
  m[[1L]] <- as.name("model.frame")
  m$... <- NULL
  mf <- eval(m, parent.frame())
  DNAME <- paste(names(mf), collapse = " by ")
  names(mf) <- NULL
  response <- attr(attr(mf, "terms"), "response")
  g <- factor(mf[[-response]])
  if (nlevels(g) != 2L)
    stop("grouping factor must have exactly 2 levels")
  DATA <- split(mf[[response]], g)
  names(DATA) <- c("x", "y")
  y <- do.call("mosesTest", c(DATA, list(...)))
  y$data.name <- DNAME
  y
  
}



#' @rdname mosesTest
#' @export
mosesTest.default <- function(x, y, extreme = NULL, ...){
  
  # example
  # x <- c(0.80, 0.83, 1.89, 1.04, 1.45, 1.38, 1.91, 1.64, 0.73, 1.46)
  # y <- c(1.15, 0.88, 0.90, 0.74, 1.21)
  # mosesTest(y, x)
  
  if(is.null(extreme)) extreme <- pmax(floor(0.05 * length(x)), 1)
  h <- extreme
  if(2*h > length(x)-2) h <- floor((length(x)-2)/2)
  
  # no alternative for the moses.test
  DNAME <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
  
  nk <- length(x)
  ne <- length(y)
  # decreasing ranks following SPSS-calculation
  R1 <- rank(-c(x, y))[1:nk]
  R1 <- sort(R1)[(h+1):(length(R1)-h)]
  
  S <- ceiling(max(R1) - min(R1) + 1)
  
  tmp <- 0
  for( i in 0 : (S - nk + 2*h)) {
    tmp <- tmp + choose(i + nk - 2*h - 2, i) * choose(ne + 2*h + 1 - i, ne - i)
  }
  
  PVAL <- (tmp / choose(nk + ne, nk))
  
  structure(list(statistic = c(S = S),
                 p.value = PVAL,
                 method = "Moses Test of Extreme Reactions",
                 alternative = "extreme values are more likely in x than in y",
                 data.name = DNAME),
            class = "htest")
  
}

