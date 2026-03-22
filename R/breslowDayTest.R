
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
#' @family topic.contingencyTests
#' @concept odds ratio
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
#' @export
breslowDayTest <- function(x, OR = NA, correct = FALSE) {
  
  if (!is.array(x) || length(dim(x)) != 3 || any(dim(x)[1:2] != 2))
    stop("x must be a 2x2xK array")
  
  K <- dim(x)[3]
  
  if (is.na(OR)) {
    a <- x[1,1,]; b <- x[1,2,]; c <- x[2,1,]; d <- x[2,2,]
    n <- a + b + c + d
    or.hat.mh <- sum(a*d/n) / sum(b*c/n)
  } else {
    or.hat.mh <- OR
  }
  
  X2.HBD <- 0
  a <- tildea <- Var.a <- numeric(K)
  
  for (j in 1:K) {
    
    mj <- rowSums(x[,,j])
    nj <- colSums(x[,,j])
    
    A <- 1 - or.hat.mh
    B <- nj[2] - mj[1] + or.hat.mh * (nj[1] + mj[1])
    C <- -mj[1] * nj[1] * or.hat.mh
    
    if (abs(A) < 1e-12) {
      tildeaj <- -C / B
    } else {
      disc <- B^2 - 4*A*C
      if (disc < 0) stop("No real solution")
      sols <- c((-B + sqrt(disc))/(2*A), (-B - sqrt(disc))/(2*A))
      tildeaj <- sols[sols > 0 & sols <= min(nj[1], mj[1])]
    }
    
    aj <- x[1,1,j]
    
    tildebj <- mj[1] - tildeaj
    tildecj <- nj[1] - tildeaj
    tildedj <- mj[2] - tildecj
    
    if (any(c(tildeaj, tildebj, tildecj, tildedj) <= 0))
      stop("Invalid expected counts")
    
    Var.aj <- 1 / (1/tildeaj + 1/tildebj + 1/tildecj + 1/tildedj)
    
    X2.HBD <- X2.HBD + (aj - tildeaj)^2 / Var.aj
    
    a[j] <- aj; tildea[j] <- tildeaj; Var.a[j] <- Var.aj
  }
  
  X2.HBDT <- X2.HBD - (sum(a) - sum(tildea))^2 / sum(Var.a)
  
  STATISTIC <- if(correct) X2.HBDT else X2.HBD
  PARAMETER <- K - 1
  PVAL <- pchisq(STATISTIC, PARAMETER, lower.tail = FALSE)
  
  structure(list(
    statistic = c("X-squared" = STATISTIC),
    parameter = c("df" = PARAMETER),
    p.value = PVAL,
    method = if(correct)
      "Breslow-Day Test (Tarone corrected)"
    else "Breslow-Day Test",
    data.name = deparse(substitute(x))
  ), class = "htest")
}



# 
# breslowDayTest <- function(x, OR = NA, correct = FALSE) {
#   
#   # Function to perform the Breslow and Day (1980) test including the
#   # corrected test by Tarone 
#   # Uses the equations in Lachin (2000),
#   # Biostatistical Methods, Wiley, p. 124-125.
#   #
#   # Programmed by Michael Hoehle <http://www.math.su.se/~hoehle>
#   # Code taken originally from a Biostatistical Methods lecture
#   # held at the Technical University of Munich in 2008.
#   #
#   # Params:
#   #  x - a 2x2xK contingency table
#   #  correct - if TRUE Tarones correction is returned
#   #
#   # Returns:
#   #  a vector with three values
#   #   statistic - Breslow and Day test statistic
#   #   pval - p value evtl. based on the Tarone test statistic
#   #               using a \chi^2(K-1) distribution
#   #
#   
#   
#   if(is.na(OR)) {
#     #Find the common OR based on Mantel-Haenszel
#     or.hat.mh <- mantelhaen.test(x)$estimate
#   } else {
#     or.hat.mh <- OR
#   }
#   
#   #Number of strata
#   K <- dim(x)[3]
#   #Value of the Statistic
#   X2.HBD <- 0
#   #Value of aj, tildeaj and Var.aj
#   a <- tildea <- Var.a <- numeric(K)
#   
#   for (j in 1:K) {
#     #Find marginals of table j
#     mj <- apply(x[,,j], MARGIN=1, sum)
#     nj <- apply(x[,,j], MARGIN=2, sum)
#     
#     #Solve for tilde(a)_j
#     coef <- c(-mj[1]*nj[1] * or.hat.mh, nj[2]-mj[1]+or.hat.mh*(nj[1]+mj[1]),
#               1-or.hat.mh)
#     sols <- Re(polyroot(coef))
#     
#     #Take the root, which fulfills 0 < tilde(a)_j <= min(n1_j, m1_j)
#     tildeaj <- sols[(0 < sols) &  (sols <= min(nj[1],mj[1]))]
#     #Observed value
#     aj <- x[1,1,j]
#     
#     #Determine other expected cell entries
#     tildebj <- mj[1] - tildeaj
#     tildecj <- nj[1] - tildeaj
#     tildedj <- mj[2] - tildecj
#     
#     #Compute \hat{\Var}(a_j | \widehat{\OR}_MH)
#     Var.aj <- (1/tildeaj + 1/tildebj + 1/tildecj + 1/tildedj)^(-1)
#     
#     #Compute contribution
#     X2.HBD <- X2.HBD + as.numeric((aj - tildeaj)^2 / Var.aj)
#     
#     #Assign found value for later computations
#     a[j] <- aj ;  tildea[j] <- tildeaj ; Var.a[j] <- Var.aj
#   }
#   
#   # Compute Tarone corrected test
#   # Add on 2015: The original equation from the 2008 lecture is incorrect
#   # as pointed out by Jean-Francois Bouzereau.
#   # X2.HBDT <-as.numeric( X2.HBD -  (sum(a) - sum(tildea))^2/sum(Var.aj) )
#   X2.HBDT <-as.numeric( X2.HBD -  (sum(a) - sum(tildea))^2/sum(Var.a) )
#   
#   DNAME <- deparse(substitute(x))
#   
#   STATISTIC <- if(correct) X2.HBDT else X2.HBD
#   PARAMETER <- K - 1
#   # Compute p-value based on the Tarone corrected test
#   PVAL <- 1 - pchisq(STATISTIC, PARAMETER)
#   METHOD <- if(correct) "Breslow-Day Test on Homogeneity of Odds Ratios (with Tarone correction)" else
#     "Breslow-Day test on Homogeneity of Odds Ratios"
#   names(STATISTIC) <- "X-squared"
#   names(PARAMETER) <- "df"
#   structure(list(statistic = STATISTIC, parameter = PARAMETER,
#                  p.value = PVAL, method = METHOD, data.name = DNAME
#   ), class = "htest")
#   
# }
# 
