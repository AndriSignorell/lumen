
#' Siegel-Tukey Test for Equality in Variability
#'
#' A nonparametric test for differences in scale (variability) between two
#' independent groups. Ranks are assigned by alternating between the extremes
#' and the center of the combined sorted sample, and a Wilcoxon-Mann-Whitney
#' statistic is then computed on these ranks.
#'
#' @description
#' The Siegel-Tukey test examines the null hypothesis that the variability
#' (scale) of \code{x} and \code{y} is equal. Rejection indicates that the
#' two groups differ in spread. The test is distribution-free and does not
#' assume normality, but it does assume equal medians under the null hypothesis
#' of equal scale. If the medians differ, the test may detect that difference
#' rather than a difference in scale; use \code{adjust.median = TRUE} to
#' remove median differences before testing.
#'
#' Ranks are assigned to the combined sorted sample in the pattern
#' 1, 4, 5, 8, 9, \ldots from the extremes inward, and 2, 3, 6, 7, \ldots
#' from the second-lowest upward. If the combined sample size is odd, the
#' median observation is dropped before ranking (it is taken from the larger
#' group when group sizes differ).
#'
#' Ties receive average ranks. The p-value is computed exactly (via
#' \code{\link{pwilcox}}) when there are no ties and both samples are smaller
#' than 50 observations; otherwise a normal approximation with tie-corrected
#' variance is used. This behaviour can be overridden with \code{exact}.
#'
#' \strong{Note:} The Siegel-Tukey test has relatively low power compared to
#' alternatives such as \code{\link{ansari.test}} or \code{\link{mood.test}},
#' and may indicate significance due to median differences rather than scale
#' differences when \code{adjust.median = FALSE}.
#'
#' @name siegelTukeyTest
#' @aliases siegelTukeyTest siegelTukeyTest.default siegelTukeyTest.formula
#'
#' @param x,y numeric vectors of data values.
#' @param adjust.median logical; if \code{TRUE}, the median of \code{x} is
#'   shifted to equal the median of \code{y} before ranking, to prevent median
#'   differences from inflating the test statistic. Default is \code{FALSE}.
#' @param alternative a character string specifying the alternative hypothesis:
#'   \code{"two.sided"} (default), \code{"greater"}, or \code{"less"}.
#'   Partial matching is allowed.
#' @param mu a single number specifying the location parameter under the null
#'   hypothesis. Default is \code{0}.
#' @param exact logical; if \code{TRUE}, an exact p-value is computed via
#'   \code{\link{pwilcox}}. Exact computation is not possible in the presence
#'   of ties; a warning is issued and the normal approximation is used instead.
#'   If \code{NA} (default), exact computation is used when both samples have
#'   fewer than 50 observations and there are no ties.
#' @param correct logical; if \code{TRUE} (default), a continuity correction
#'   is applied in the normal approximation. Ignored when \code{exact = TRUE}
#'   or when ties are present (continuity correction is not appropriate with
#'   tie-corrected variance).
#' @param formula a formula of the form \code{response ~ group}, where
#'   \code{response} is a numeric vector and \code{group} a factor or vector
#'   with exactly two levels.
#' @param data an optional data frame (or matrix, coerced to data frame)
#'   containing the variables in \code{formula}. If not supplied, variables
#'   are taken from \code{environment(formula)}.
#' @param subset an optional vector specifying a subset of observations to use.
#' @param na.action a function indicating how to handle \code{NA}s in the
#'   formula interface. Defaults to \code{na.pass}; \code{NA}s in \code{x} or
#'   \code{y} are silently dropped in the default method.
#' @param \dots further arguments passed to or from methods.
#'
#' @return An object of class \code{"htest"} with the following components:
#' \describe{
#'   \item{\code{statistic}}{the Wilcoxon rank-sum statistic \eqn{W} computed
#'     on the Siegel-Tukey ranks.}
#'   \item{\code{p.value}}{the p-value of the test.}
#'   \item{\code{null.value}}{the location parameter \code{mu} under the null
#'     hypothesis.}
#'   \item{\code{alternative}}{a character string describing the alternative
#'     hypothesis.}
#'   \item{\code{method}}{a character string identifying the test.}
#'   \item{\code{data.name}}{a character string giving the names of the data.}
#'   \item{\code{exact}}{logical indicating whether the exact distribution was
#'     used.}
#'   \item{\code{ties}}{logical indicating whether ties were present in the
#'     Siegel-Tukey ranks.}
#' }
#'
#' @seealso \code{\link{ansari.test}}, \code{\link{mood.test}},
#'   \code{\link{wilcox.test}}, \code{\link{leveneTest}}
#'
#' @references
#' Siegel, S. and Tukey, J. W. (1960): A nonparametric sum of ranks procedure
#' for relative spread in unpaired samples. \emph{Journal of the American
#' Statistical Association}, \bold{55}(291), 429--445.
#'
#' Sheskin, D. J. (2004): \emph{Handbook of Parametric and Nonparametric
#' Statistical Procedures}, 3rd ed. Chapman & Hall/CRC, Boca Raton, FL.
#'
#' @note
#' Originally based on a blog post by Tal Galili:
#' \url{https://www.r-statistics.com/2010/02/siegel-tukey-a-non-parametric-test-for-equality-in-variability-r-code/}
#'
#' @examples
#' # Duller, S. 183
#' x <- c(12, 13, 29, 30)
#' y <- c(15, 17, 18, 24, 25, 26)
#' siegelTukeyTest(x, y)
#' siegelTukeyTest(x, y, alternative = "greater")
#'
#' # Duller, S. 323
#' old <- c(870, 930, 935, 1045, 1050, 1052, 1055)
#' new <- c(932, 970, 980, 1001, 1009, 1030, 1032, 1040, 1046)
#' siegelTukeyTest(old, new, alternative = "greater")
#' # recommended alternatives:
#' mood.test(old, new, alternative = "greater")
#' ansari.test(old, new, alternative = "greater")
#'
#' # Bortz, S. 250 -- formula interface
#' x  <- c(26.3, 26.5, 26.8, 27.0, 27.0, 27.2, 27.3,
#'         27.3, 27.4, 27.5, 27.6, 27.8, 27.9)
#' id <- c(2, 2, 2, 1, 2, 2, 1, 2, 2, 1, 1, 1, 2) - 1
#' siegelTukeyTest(x ~ id)
#'
#' # Sachs (2007), S. 314
#' A <- c(10.1, 7.3, 12.6, 2.4, 6.1, 8.5, 8.8, 9.4, 10.1, 9.8)
#' B <- c(15.3, 3.6, 16.5, 2.9, 3.3, 4.2, 4.9, 7.3, 11.7, 13.1)
#' siegelTukeyTest(A, B)
#'
#' # equal medians, different spread
#' x <- c(4, 4, 5, 5, 6, 6)
#' y <- c(0, 0, 1, 9, 10, 10)
#' siegelTukeyTest(x, y)
#'
#' # unequal group sizes
#' x <- c(4, 4, 5, 5, 6, 6)
#' y <- c(0, 0, 1, 9, 10)
#' siegelTukeyTest(x, y)
#'
#' # median adjustment
#' x <- c(177, 200, 227, 230, 232, 268, 272, 297)
#' y <- c(47, 105, 126, 142, 158, 172, 197, 220, 225, 230, 262, 270)
#' siegelTukeyTest(x, y, adjust.median = TRUE)
#'
#' # floating-point robustness (previously affected by merge precision bug)
#' y2 <- c(-1, 2, 2.1, 3)
#' x2 <- c(-5, -9, 13, 12, 90, 100)
#' siegelTukeyTest(x2, y2, adjust.median = TRUE)  
#' # p ~ 0.1143
#' 

#' @rdname siegelTukeyTest
#' @family test.variance
#' @concept hypothesis-testing
#' @concept nonparametric
#'
#'


#' @export
siegelTukeyTest <- function (x, ...)  UseMethod("siegelTukeyTest")

# compare:  jmuOutlier::siegel.test()



#' @rdname siegelTukeyTest
#' @export
siegelTukeyTest.formula <- function(formula, data, subset,
                                    na.action = na.pass, ...) {
  
  if (missing(formula) || length(formula) != 3L)
    stop("'formula' missing or incorrect")
  
  args <- list(
    formula   = formula,
    na.action = na.action,
    allowed   = "two.sample.independent"
  )
  
  if (!missing(data))
    args$data <- data
  
  if (!missing(subset))
    args$subset <- substitute(subset)
  
  d <- do.call(bedrock::resolveFormula, args)
  
  siegelTukeyTest.default(
    x = d$x,
    y = d$y,
    ...
  )
}



#' @rdname siegelTukeyTest
#' @export
siegelTukeyTest.default <- function(x, y, alternative = c("two.sided", "less", "greater"),
                                    mu = 0, adjust.median = FALSE, exact = NA,
                                    correct = TRUE, ...) {
  
  alternative <- match.arg(alternative)
  
  # validate logical arguments
  if (!is.logical(adjust.median) || length(adjust.median) != 1L || is.na(adjust.median))
    stop("'adjust.median' must be TRUE or FALSE")
  
  if (!is.logical(correct) || length(correct) != 1L || is.na(correct))
    stop("'correct' must be TRUE or FALSE")
  
  if (!is.logical(exact) || length(exact) != 1L)
    stop("'exact' must be TRUE, FALSE, or NA")
  
  if (!missing(mu) && ((length(mu) > 1L) || !is.finite(mu)))
    stop("'mu' must be a single number")
  
  if (missing(y))
    stop("'y' is missing")
  
  if (!is.numeric(x) || !is.numeric(y))
    stop("'x' and 'y' must be numeric")
  
  # remove non-finite values
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  
  if (length(x) < 2L || length(y) < 2L)
    stop("'x' and 'y' must contain at least two observations each")
  
  dname <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
  
  if (adjust.median)
    x <- x - (median(x) - median(y))
  
  xx <- c(x, y)
  id <- c(rep(0, length(x)), rep(1, length(y)))
  
  strank <- .siegelTukeyRank(xx, g = id)
  
  # ties defined by duplicated raw values, not by coincidentally equal mean ranks
  TIES <- anyDuplicated(strank$sort.x) > 0
  
  # honour explicit exact argument, auto-detect only when NA
  if (is.na(exact))
    exact <- (length(x) < 50L) && (length(y) < 50L) && !TIES
  
  if (exact && TIES) {
    warning("cannot compute exact p-value with ties, using normal approximation")
    exact <- FALSE
  }
  
  ranks0 <- strank$unique.ranks[strank$sort.id == 0]
  ranks1 <- strank$unique.ranks[strank$sort.id == 1]
  
  m <- length(ranks0)
  n <- length(ranks1)
  N <- m + n
  
  W <- sum(ranks1)
  U <- W - n * (n + 1) / 2
  
  if (exact) {
    # exact Wilcoxon distribution applied to Siegel-Tukey ranks (valid without ties)
    p.value <- switch(alternative,
                      two.sided = {
                        p <- if (U > (m * n / 2))
                          pwilcox(U - 1, m, n, lower.tail = FALSE)
                        else
                          pwilcox(U, m, n)
                        min(2 * p, 1)
                      },
                      greater = pwilcox(U - 1, m, n, lower.tail = FALSE),
                      less    = pwilcox(U, m, n)
    )
    
  } else {
    # normal approximation with tie-corrected variance
    tie_table     <- table(strank$unique.ranks)
    tie_correction <- sum(tie_table^3 - tie_table) / (N * (N - 1) * (N + 1))
    
    E <- m * n / 2
    V <- m * n * (N + 1) / 12 * (1 - tie_correction)
    
    z <- U - E
    
    if (correct) {
      CORRECTION <- switch(alternative,
                           two.sided = sign(z) * 0.5,
                           greater   =  0.5,
                           less      = -0.5
      )
      z <- z - CORRECTION
    }
    
    z <- z / sqrt(V)
    
    p.value <- switch(alternative,
                      less      = pnorm(z),
                      greater   = pnorm(z, lower.tail = FALSE),
                      two.sided = 2 * min(pnorm(z), pnorm(z, lower.tail = FALSE))
    )
  }
  
  # mu is accepted for interface compatibility but not used in the test statistic;
  # the ST test has no natural location shift parameter
  structure(
    list(
      statistic   = c(W = W),
      parameter   = NULL,
      p.value     = p.value,
      null.value  = c(mu = mu),
      alternative = alternative,
      method      = "Siegel-Tukey test for scale differences",
      data.name   = dname,
      exact       = exact,
      ties        = TIES
    ),
    class = "htest"
  )
}


# == internal helper functions ==============================================

.siegelTukeyRank <- function(x, g, dropMedian = TRUE) {
  
  ord.x <- order(x, g)
  sort.x <- x[ord.x]
  sort.id <- g[ord.x]
  
  n <- length(x)
  if(dropMedian){
    if(n %% 2 > 0) {
      fm <- which(sort.x == median(sort.x))[1]
      sort.x <- sort.x[-fm]
      sort.id <- sort.id[-fm]
      n <- n - 1
    }
  }
  
  base1 <- c(1, 4)
  iterator1 <- matrix(seq(from = 1, to = n, by = 4)) - 1
  rank1 <- apply(iterator1, 1, function(x) x + base1)
  
  iterator2 <- matrix(seq(from = 2, to = n, by = 4))
  base2 <- c(0, 1)
  rank2 <- apply(iterator2, 1, function(x) x + base2)
  
  if (length(rank1) == length(rank2)) {
    rank <- c(rank1[1:floor(n/2)], rev(rank2[1:ceiling(n/2)]))
  } else {
    rank <- c(rank1[1:ceiling(n/2)], rev(rank2[1:floor(n/2)]))
  }
  
  # NEW: Match using a positional index instead of a merge, 
  #      avoiding floating-point comparisons
  grp <- match(sort.x, unique(sort.x))
  unique.ranks <- tapply(rank, grp, mean)
  avg.ranks <- unique.ranks[grp]
  
  data.frame(
    sort.x       = sort.x,
    sort.id      = sort.id,
    unique.ranks = as.numeric(avg.ranks),
    raw.ranks    = rank
  )
  
}
