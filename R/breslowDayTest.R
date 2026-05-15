
#' Breslow-Day Test for Homogeneity of the Odds Ratios 
#' 
#' A test for homogeneity of odds ratios across several 2×2 contingency 
#' tables (strata), commonly used to verify whether a confounding effect 
#' is constant across subgroups, as required by the Mantel-Haenszel method.

#' Calculates the Breslow-Day test of homogeneity for a \eqn{2 \times 2 \times
#' k}{2 x 2 x k} table, in order to investigate if all \eqn{k} strata have the
#' same OR. If OR is not given, the Mantel-Haenszel estimate is used. 
#' 
#' For the Breslow-Day test to be valid, the sample size should be relatively
#' large in each stratum, and at least 80\% of the expected cell counts should
#' be greater than 5. Note that this is a stricter sample size requirement than
#' the requirement for the Cochran-Mantel-Haenszel test for tables, in that
#' each stratum sample size (not just the overall sample size) must be
#' relatively large. Even when the Breslow-Day test is valid, it might not be
#' very powerful against certain alternatives, as discussed in Breslow and Day
#' (1980).
#' 
#' Alternatively, it might be better to cast the entire inference problem into
#' the setting of a logistic regression model. Here, the underlying question of
#' the Breslow-Day test can be answered by investigating whether an interaction
#' term with the strata variable is necessary (e.g. using a likelihood ratio
#' test using the \code{anova} function).
#' 
#' @name breslowDayTest
#' @param x a \eqn{2 \times 2 \times k}{2 x 2 x k} table. 
#' @param OR the odds ratio to be tested against. If left undefined (default)
#' the Mantel-Haenszel estimate will be used. 
#' @param correct If TRUE, the Breslow-Day test with Tarone's adjustment is
#' computed, which subtracts an adjustment factor to make the resulting
#' statistic asymptotically chi-square. 
#' 
#' @note
#' Based on code by Michael Hoehle.
#' 
#' @seealso \code{\link{woolfTest}} 
#' @references Breslow, N. E., N. E. Day (1980) The Analysis of Case-Control
#' Studies \emph{Statistical Methods in Cancer Research: Vol. 1}. Lyon, France,
#' IARC Scientific Publications.
#' 
#' Tarone, R.E. (1985) On heterogeneity tests based on efficient scores,
#' \emph{Biometrika}, 72, pp. 91-95.
#' 
#' Jones, M. P., O'Gorman, T. W., Lemka, J. H., and Woolson, R. F. (1989) A
#' Monte Carlo Investigation of Homogeneity Tests of the Odds Ratio Under
#' Various Sample Size Configurations \emph{Biometrics}, 45, 171-181 \cr
#' 
#' Breslow, N. E. (1996) Statistics in Epidemiology: The Case-Control Study
#' \emph{Journal of the American Statistical Association}, 91, 14-26.
#' 
#' @examples
#' 
#' migraine <- xtabs(freq ~ .,
#'             cbind(expand.grid(treatment=c("active", "placebo"),
#'                               response =c("better", "same"),
#'                               gender   =c("female", "male")),
#'                   freq=c(16, 5, 11, 20, 12, 7, 16, 19))
#'             )
#' 
#' # get rid of gender
#' tab <- xtabs(Freq ~ treatment + response, migraine)
#' tab
#' 
#' # only the women
#' female <- migraine[,, 1]
#' female
#' 
#' # .. and the men
#' male <- migraine[,, 2]
#' male
#' 
#' breslowDayTest(migraine)
#' breslowDayTest(migraine, correct = TRUE)
#' 
#' 
#' salary <- array(
#'       c(38, 12, 102, 141, 12, 9, 136, 383),
#'       dim=c(2, 2, 2),
#'       dimnames=list(exposure=c("exposed", "not"),
#'                     disease =c("case", "control"),
#'                     salary  =c("<1000", ">=1000"))
#'                     )
#' 
#' # common odds ratio = 4.028269
#' breslowDayTest(salary, OR = 4.02)
#' 

#' @rdname breslowDayTest
#' @family test.contingency
#' @concept hypothesis-testing
#' @concept table-manipulation
#'
#'


#' @export
breslowDayTest <- function(x, OR = NA, correct = FALSE) {
  
  ## -------------------------------------------------------------------
  ## Input validation
  ## -------------------------------------------------------------------
  if (!is.array(x) || length(dim(x)) != 3L || any(dim(x)[1:2] != 2L))
    stop("'x' must be a 2x2xK array")
  
  if (any(x < 0, na.rm = TRUE) || any(!is.finite(x)))
    stop("all entries of 'x' must be nonnegative and finite")
  
  if (any(x != round(x)))
    warning("'x' contains non-integer counts", call. = FALSE)
  
  correct <- as.logical(correct)
  if (length(correct) != 1L || is.na(correct))
    stop("'correct' must be TRUE or FALSE")
  
  if (!is.na(OR)) {
    if (!is.numeric(OR) || length(OR) != 1L || !is.finite(OR) || OR <= 0)
      stop("'OR' must be a positive finite number")
  }
  
  K <- dim(x)[3L]
  
  ## -------------------------------------------------------------------
  ## Mantel-Haenszel estimate of common OR
  ## -------------------------------------------------------------------
  if (is.na(OR)) {
    a <- x[1,1,]; b <- x[1,2,]; c <- x[2,1,]; d <- x[2,2,]
    n <- a + b + c + d
    denom <- sum(b * c / n)
    if (denom == 0)
      stop("Mantel-Haenszel denominator is zero; cannot estimate common OR")
    or.hat.mh <- sum(a * d / n) / denom
  } else {
    or.hat.mh <- OR
  }
  
  ## -------------------------------------------------------------------
  ## Per-stratum computation
  ## -------------------------------------------------------------------
  X2.HBD <- 0
  a <- tildea <- Var.a <- numeric(K)
  
  for (j in seq_len(K)) {
    
    mj <- rowSums(x[,,j])
    nj <- colSums(x[,,j])
    
    if (any(mj == 0) || any(nj == 0))
      stop("stratum ", j, " has zero marginal totals; ",
           "the test is not defined", call. = FALSE)
    
    A <- 1 - or.hat.mh
    B <- nj[2] - mj[1] + or.hat.mh * (nj[1] + mj[1])
    C <- -mj[1] * nj[1] * or.hat.mh
    
    if (abs(A) < 1e-12) {
      tildeaj <- -C / B
    } else {
      disc <- B^2 - 4 * A * C
      ## Guard against small negative disc from floating point
      if (disc < -sqrt(.Machine$double.eps))
        stop("no real solution for expected count in stratum ", j,
             call. = FALSE)
      disc    <- max(disc, 0)
      sols    <- c((-B + sqrt(disc)) / (2*A), (-B - sqrt(disc)) / (2*A))
      tildeaj <- sols[sols > 0 & sols <= min(nj[1], mj[1])]
    }
    
    if (length(tildeaj) == 0L)
      stop("no valid solution for expected count in stratum ", j,
           call. = FALSE)
    if (length(tildeaj) > 1L)
      stop("ambiguous solution for expected count in stratum ", j,
           "; this should not occur with valid data", call. = FALSE)
    
    aj      <- x[1,1,j]
    tildebj <- mj[1] - tildeaj
    tildecj <- nj[1] - tildeaj
    tildedj <- mj[2] - tildecj
    
    if (any(c(tildeaj, tildebj, tildecj, tildedj) <= 0))
      warning("non-positive expected counts in stratum ", j,
              "; results may be unreliable", call. = FALSE)
    
    Var.aj  <- 1 / (1/tildeaj + 1/tildebj + 1/tildecj + 1/tildedj)
    X2.HBD  <- X2.HBD + (aj - tildeaj)^2 / Var.aj
    
    a[j]      <- aj
    tildea[j] <- tildeaj
    Var.a[j]  <- Var.aj
  }
  
  ## -------------------------------------------------------------------
  ## Tarone correction
  ## -------------------------------------------------------------------
  if (correct) {
    if (sum(Var.a) <= 0)
      stop("sum of variances is zero; Tarone correction is undefined",
           call. = FALSE)
    X2.HBDT <- X2.HBD - (sum(a) - sum(tildea))^2 / sum(Var.a)
  }
  
  STATISTIC <- unname(if (correct) X2.HBDT else X2.HBD)
  PARAMETER <- K - 1L
  
  structure(
    list(
      statistic = c("X-squared" = STATISTIC),
      parameter = c(df = PARAMETER),
      p.value   = pchisq(STATISTIC, PARAMETER, lower.tail = FALSE),
      method    = if (correct)
        "Breslow-Day Test (Tarone corrected)"
      else
        "Breslow-Day Test",
      data.name = deparse(substitute(x)),
      n         = sum(x)
    ),
    class = "htest"
  )
}

