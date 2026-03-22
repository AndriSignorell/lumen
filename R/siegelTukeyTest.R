

#' Siegel-Tukey Test For Equality In Variability
#' 
#' A nonparametric test for differences in scale (variability) between two 
#' independent groups, based on a specific rank assignment that alternates 
#' between the extremes and the center of the combined sample.
#' 
#' Non-parametric Siegel-Tukey test for equality in variability. The null
#' hypothesis is that the variability of x is equal between two groups. A
#' rejection of the null hypothesis indicates that variability differs between
#' the two groups. 
#' 
#' \strong{Note:} \verb{   } The Siegel-Tukey test has relatively low power 
#' and may, under certain
#' conditions, indicate significance due to differences in medians rather than
#' differences in variabilities (consider using the argument
#' \code{adjust.median}). 
#' 
#' @name siegelTukeyTest
#' @aliases siegelTukeyTest siegelTukeyTest.default siegelTukeyTest.formula
#' 
#' @param x,y numeric vector of data values. Non-finite (e.g. infinite or
#' missing) values will be omitted.
#' @param adjust.median Should between-group differences in medians be leveled
#' before performing the test? In certain cases, the Siegel-Tukey test is
#' susceptible to median differences and may indicate significant differences
#' in variability that, in reality, stem from differences in medians. Default
#' is \code{FALSE}.
#' @param alternative a character string specifying the alternative hypothesis,
#' must be one of \code{"two.sided"} (default), \code{"greater"} or
#' \code{"less"}. You can specify just the initial letter.
#' @param mu a number specifying an optional parameter used to form the null
#' hypothesis. See Details.
#' @param exact a logical indicating whether an exact p-value should be
#' computed. This is passed directly to \code{\link{wilcox.test}}.
#' @param correct a logical indicating whether to apply continuity correction
#' in the normal approximation for the p-value.
#' @param conf.int a logical indicating whether a confidence interval should be
#' computed.
#' @param conf.level confidence level of the interval.
#' @param formula a formula of the form \code{lhs ~ rhs} where \code{lhs} gives
#' the data values and rhs the corresponding groups.
#' @param data an optional matrix or data frame (or similar: see
#' \code{\link{model.frame}}) containing the variables in the formula
#' \code{formula}.  By default the variables are taken from
#' \code{environment(formula)}.
#' @param subset an optional vector specifying a subset of observations to be
#' used.
#' @param na.action a function which indicates what should happen when the data
#' contain NAs. Defaults to \code{getOption("na.action")}.
#' @param \dots further arguments to be passed to or from methods.
#' @return A list of class \code{htest}, containing the following components:
#' \item{statistic}{ Siegel-Tukey test (Wilcoxon test on tie-adjusted
#' Siegel-Tukey ranks, after the median adjustment if specified).}
#' \item{p.value}{ the p-value for the test} \item{null.value}{is the value of
#' the median specified by the null hypothesis. This equals the input argument
#' \code{mu}. } \item{alternative}{a character string describing the
#' alternative hypothesis.} \item{method}{ the type of test applied}
#' \item{data.name}{a character string giving the names of the data.}
#' 
#' @note
#' Based on a blog post by Daniel Malter, Tal Galili <tal.galili@@gmail.com>
#' 
#' \href{https://www.r-statistics.com/2010/02/siegel-tukey-a-non-parametric-test-for-equality-in-variability-r-code/}{www.r-statistics.com/2010/02/}
#' 
#' @seealso \code{\link{mood.test}}, \code{\link{ansari.test}},
#' \code{\link{wilcox.test}}, \code{\link{leveneTest}}
#' 
#' @references 
#' Sheskin, D. J. (2004): \emph{Handbook of parametric and nonparametric
#' statistical procedures} 3rd edition. Chapman and Hall/CRC. Boca Raton, FL.
#' 
#' Siegel, S., Tukey, J. W. (1960): A nonparametric sum of ranks
#' procedure for relative spread in unpaired samples. \emph{Journal of the
#' American Statistical Association}.
#' 
#' @family topic.dispersionTests
#' @family topic.nonparametricTests
#' @concept scale test
#' 
#' @examples
#' 
#' # Duller, S. 183
#' x <- c(12, 13, 29, 30)
#' y <- c(15, 17, 18, 24, 25, 26)
#' siegelTukeyTest(x, y)
#' siegelTukeyTest(x, y, alternative="greater")
#' 
#' # Duller, S. 323
#' old <- c(870,930,935,1045,1050,1052,1055)
#' new <- c(932,970,980,1001,1009,1030,1032,1040,1046)
#' siegelTukeyTest(old, new, alternative = "greater")
#' # compare to the recommended alternatives
#' mood.test(old, new, alternative="greater")
#' ansari.test(old, new, alternative="greater")
#' 
#' # Bortz, S. 250
#' x <- c(26.3,26.5,26.8,27.0,27.0,27.2,27.3,27.3,27.4,27.5,27.6,27.8,27.9)
#' id <- c(2,2,2,1,2,2,1,2,2,1,1,1,2)-1
#' siegelTukeyTest(x ~ id)
#' 
#' 
#' # Sachs, Angewandte Statistik, 12. Auflage, 2007, S. 314
#' A <- c(10.1,7.3,12.6,2.4,6.1,8.5,8.8,9.4,10.1,9.8)
#' B <- c(15.3,3.6,16.5,2.9,3.3,4.2,4.9,7.3,11.7,13.1)
#' siegelTukeyTest(A, B)
#' 
#' 
#' 
#' ### 1
#' x <- c(4,4,5,5,6,6)
#' y <- c(0,0,1,9,10,10)
#' siegelTukeyTest(x, y)
#' 
#' ### 2
#' # example for a non equal number of cases:
#' x <- c(4,4,5,5,6,6)
#' y <- c(0,0,1,9,10)
#' siegelTukeyTest(x, y)
#' 
#' ### 3
#' x <- c(33, 62, 84, 85, 88, 93, 97, 4, 16, 48, 51, 66, 98)
#' id <- c(0,0,0,0,0,0,0,1,1,1,1,1,1)
#' siegelTukeyTest(x ~ id)
#' 
#' ### 4
#' x <- c(177,200,227,230,232,268,272,297,47,105,126,142,158,172,197,220,225,230,262,270)
#' id <- c(rep(0,8),rep(1,12))
#' siegelTukeyTest(x ~ id, adjust.median=TRUE)
#' 
#' ### 5
#' x <- c(33,62,84,85,88,93,97)
#' y <- c(4,16,48,51,66,98)
#' siegelTukeyTest(x, y)
#' 
#' ### 6
#' x <- c(0,0,1,4,4,5,5,6,6,9,10,10)
#' id <- c(0,0,0,1,1,1,1,1,1,0,0,0)
#' siegelTukeyTest(x ~ id)
#' 
#' ### 7
#' x <- c(85,106,96, 105, 104, 108, 86)
#' id <- c(0,0,1,1,1,1,1)
#' siegelTukeyTest(x ~ id)
#' 


#' @rdname siegelTukeyTest
#' @export
siegelTukeyTest <- function (x, ...)  UseMethod("siegelTukeyTest")

# compare:  jmuOutlier::siegel.test()


#' @rdname siegelTukeyTest
#' @export
siegelTukeyTest.formula <- local({
  
  # super elegant formula implementation
  # in fact we need nothing other, than is already implemented 

  tf <- getS3method("wilcox.test", "formula")
  
  new_body <- .replace_text_calls(body(tf), old="wilcox.test", 
                                  new="siegelTukeyTest")
  
  new_fun <- tf
  body(new_fun) <- new_body
  
  new_fun
  
})





#' @rdname siegelTukeyTest
#' @export
siegelTukeyTest.default <- function(x, y, adjust.median = FALSE,
                                    alternative = c("two.sided","less","greater"), mu = 0,
                                    exact = NULL, correct = TRUE, 
                                    conf.int = FALSE, conf.level = 0.95, ...) {
  ###### published on:
  #   http://www.r-statistics.com/2010/02/siegel-tukey-a-non-parametric-test-for-equality-in-variability-r-code/
  #   Main author of the function:  Daniel Malter
  
  # Doku: http://www.crcnetbase.com/doi/abs/10.1201/9781420036268.ch14
  
  
  if (!missing(mu) && ((length(mu) > 1L) || !is.finite(mu)))
    stop("'mu' must be a single number")
  
  if (conf.int) {
    if (!((length(conf.level) == 1L) && is.finite(conf.level) &&
          (conf.level > 0) && (conf.level < 1)))
      stop("'conf.level' must be a single number between 0 and 1")
  }
  
  if (!is.numeric(x))
    stop("'x' must be numeric")
  
  if (!is.null(y)) {
    if (!is.numeric(y))
      stop("'y' must be numeric")
    DNAME <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
    x <- x[is.finite(x)]
    y <- y[is.finite(y)]
  }
  else {
    DNAME <- deparse(substitute(x))
    x <- x[is.finite(x)]
  }
  
  # adjusting median
  if (adjust.median) {
    x <- x - median(x)
    y <- y - median(y)
  }
  
  # the larger group comes first
  if( length(x) > length(y) ){
    xx <- c(x, y)
    id <- c(rep(0, length(x)), rep(1, length(y)))
  } else {
    xx <- c(y,x)
    id <- c(rep(0, length(y)), rep(1, length(x)))
  }
  
  strank <- .siegelTukeyRank(xx, g = id)
  ranks0 <- strank$unique.ranks[strank$sort.id == 0]
  ranks1 <- strank$unique.ranks[strank$sort.id == 1]
  
  RVAL <- wilcox.test(ranks0, ranks1, alternative = alternative,
                      mu = mu, paired = FALSE, exact = exact, correct = correct,
                      conf.int = conf.int, conf.level = conf.level)
  
  RVAL$statistic <- sum(ranks1)
  names(RVAL$statistic)  <- "ST"
  RVAL$data.name <- DNAME
  RVAL <- c(RVAL, list(stranks = strank, MeanRanks = c(mean(ranks0), mean(ranks1))))
  RVAL$method <- "Siegel-Tukey-test for equal variability"
  RVAL$null.value <- 1
  names(RVAL$null.value) <- "ratio of scales"
  class(RVAL) <- "htest"
  return(RVAL)
  
  if(suppressWarnings(wilcox.test(x,y)$p.value) < 0.05) 
    warning("siegelTukeyTest: wilcox.test(x, y) is significant! Consider setting adjust.median = TRUE." )
  
}




# == internal helper functions ==============================================

.siegelTukeyRank <- function(x, g, drop.median = TRUE) {
  
  # they do not drop the median in:
  # http://en.wikipedia.org/wiki/Siegel%E2%80%93Tukey_test
  # A <- c(33,62,84,85,88,93,97); B <- c(4,16,48,51,66,98)
  # this is wrong there, as the author did not leave the median out
  
  ord.x <- order(x, g)
  sort.x <- x[ord.x]
  sort.id <- g[ord.x]
  
  n <- length(x)
  if(drop.median){
    if(n %% 2 > 0) {
      # gonna have to drop the (first) median value
      # as we sorted by the groupsize, this will be the one out of the bigger group (if existing)
      fm <- which( sort.x == median(sort.x))[1]
      sort.x <- sort.x[-fm]
      sort.id <- sort.id[-fm]
      n <- n-1
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
  
  unique.ranks <- tapply(rank, sort.x, mean)
  
  # changed by 0.99.60 in consequence of https://github.com/AndriSignorell/DescTools/issues/171
  # names might not be exact...
  # unique.x <- as.numeric(as.character(names(unique.ranks)))
  unique.x <- unique(sort.x)
  
  ST.matrix <- merge(
    data.frame(sort.x, sort.id),          # this are the original values in x-order
    data.frame(unique.x, unique.ranks),   # this is the rank.matrix
    by.x = "sort.x", by.y = "unique.x")
  
  return(ST.matrix)
}



