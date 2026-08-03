
#' Breslow-Day Test  for Comparing Odds Ratios Across Strata
#'
#' A test for homogeneity of odds ratios across several 2x2 contingency
#' tables (strata), commonly used to verify whether a confounding effect
#' is constant across subgroups, as required by the Mantel-Haenszel method.
#'
#' Calculates the Breslow-Day test of homogeneity for a
#' \eqn{2 \times 2 \times k}{2 x 2 x k} table, in order to investigate if
#' all \eqn{k} strata have the same OR. If \code{OR} is not given, the
#' Mantel-Haenszel estimate is used.
#'
#' For the Breslow-Day test to be valid, the sample size should be
#' relatively large in each stratum, and at least 80\% of the expected cell
#' counts should be greater than 5. Note that this is a stricter sample
#' size requirement than the requirement for the Cochran-Mantel-Haenszel
#' test for tables, in that each stratum sample size (not just the overall
#' sample size) must be relatively large. Even when the Breslow-Day test is
#' valid, it might not be very powerful against certain alternatives, as
#' discussed in Breslow and Day (1980).
#'
#' The statistic is referred to a chi-squared distribution with \eqn{k-1}
#' degrees of freedom; this also applies when a prespecified \code{OR} is
#' supplied. Note that Tarone's adjustment is derived for the
#' Mantel-Haenszel estimate; a warning is issued if \code{correct = TRUE}
#' is combined with a user-supplied \code{OR}.
#'
#' Alternatively, it might be better to cast the entire inference problem
#' into the setting of a logistic regression model. Here, the underlying
#' question of the Breslow-Day test can be answered by investigating
#' whether an interaction term with the strata variable is necessary (e.g.
#' using a likelihood ratio test using the \code{anova} function).
#'
#' @param x a \eqn{2 \times 2 \times k}{2 x 2 x k} table.
#' @param OR the odds ratio to be tested against. If left undefined
#' (default) the Mantel-Haenszel estimate will be used.
#' @param correct logical, if \code{TRUE} the Breslow-Day test with
#' Tarone's adjustment is computed, which subtracts an adjustment factor to
#' make the resulting statistic asymptotically chi-squared.
#' @return A list with class \code{"htest"} containing the following
#' components:
#'   \item{\code{statistic}}{the value of the chi-squared test statistic.}
#'   \item{\code{parameter}}{the degrees of freedom of the approximate
#'     chi-squared distribution of the test statistic (\eqn{k-1}).}
#'   \item{\code{p.value}}{the p-value of the test.}
#'   \item{\code{method}}{a character string indicating the test performed.}
#'   \item{\code{data.name}}{a character string giving the name of the data.}
#'   \item{\code{n}}{the total number of observations (not shown on screen).}
#'
#' @note
#' Based on code by Michael Hoehle, adapted to conform to package
#' standards.
#'
#' @seealso \code{\link{mantelhaen.test}}
#'
#' @references Breslow, N. E. and Day, N. E. (1980) The Analysis of
#' Case-Control Studies. \emph{Statistical Methods in Cancer Research:
#' Vol. 1}. Lyon, France, IARC Scientific Publications.
#'
#' Tarone, R.E. (1985) On heterogeneity tests based on efficient scores,
#' \emph{Biometrika}, 72, pp. 91-95.
#'
#' Jones, M. P., O'Gorman, T. W., Lemka, J. H., and Woolson, R. F. (1989)
#' A Monte Carlo Investigation of Homogeneity Tests of the Odds Ratio
#' Under Various Sample Size Configurations. \emph{Biometrics}, 45,
#' 171-181.
#'
#' Breslow, N. E. (1996) Statistics in Epidemiology: The Case-Control
#' Study. \emph{Journal of the American Statistical Association}, 91,
#' 14-26.
#'
#' @examples
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
#' @family test.categorical
#' @concept categorical-test
#' @concept homogeneity
#'
#' @export
breslowDayTest <- function(x, OR = NULL, correct = FALSE) {

  DNAME <- deparse1(substitute(x))

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

  # OR = NULL (or NA for backward compatibility) requests the MH estimate
  if (!is.null(OR) && length(OR) == 1L && is.na(OR))
    OR <- NULL
  if (!is.null(OR)) {
    if (!is.numeric(OR) || length(OR) != 1L || !is.finite(OR) || OR <= 0)
      stop("'OR' must be a positive finite number")
    if (correct)
      warning("Tarone's adjustment is derived for the Mantel-Haenszel ",
              "estimate; interpret the corrected statistic with caution ",
              "when 'OR' is supplied", call. = FALSE)
  }

  K <- dim(x)[3L]

  ## -------------------------------------------------------------------
  ## Mantel-Haenszel estimate of common OR
  ## -------------------------------------------------------------------
  if (is.null(OR)) {
    n <- apply(x, 3L, sum)
    denom <- sum(x[1, 2, ] * x[2, 1, ] / n)
    if (denom == 0)
      stop("Mantel-Haenszel denominator is zero; cannot estimate common OR")
    or.hat.mh <- sum(x[1, 1, ] * x[2, 2, ] / n) / denom
  } else {
    or.hat.mh <- OR
  }

  ## -------------------------------------------------------------------
  ## Per-stratum computation
  ## -------------------------------------------------------------------
  X2.HBD <- 0
  a <- tildea <- Var.a <- numeric(K)

  for (j in seq_len(K)) {

    mj <- rowSums(x[, , j])
    nj <- colSums(x[, , j])
    Nj <- sum(mj)

    if (any(mj == 0) || any(nj == 0))
      stop("stratum ", j, " has zero marginal totals; ",
           "the test is not defined", call. = FALSE)

    A <- 1 - or.hat.mh
    B <- nj[2] - mj[1] + or.hat.mh * (nj[1] + mj[1])
    C <- -mj[1] * nj[1] * or.hat.mh

    # admissible range for the expected count of cell (1,1)
    lo <- max(0, mj[1] + nj[1] - Nj)
    hi <- min(nj[1], mj[1])

    if (abs(A) < 1e-12) {
      tildeaj <- -C / B
    } else {
      disc <- B^2 - 4 * A * C
      ## Guard against small negative disc from floating point
      if (disc < -sqrt(.Machine$double.eps))
        stop("no real solution for expected count in stratum ", j,
             call. = FALSE)
      disc    <- max(disc, 0)
      sols    <- c((-B + sqrt(disc)) / (2 * A),
                   (-B - sqrt(disc)) / (2 * A))
      tildeaj <- sols[sols > lo & sols < hi]
    }

    if (length(tildeaj) == 0L)
      stop("no valid solution for expected count in stratum ", j,
           call. = FALSE)
    if (length(tildeaj) > 1L)
      stop("ambiguous solution for expected count in stratum ", j,
           "; this should not occur with valid data", call. = FALSE)

    aj      <- x[1, 1, j]
    tildebj <- mj[1] - tildeaj
    tildecj <- nj[1] - tildeaj
    tildedj <- mj[2] - tildecj

    Var.aj  <- 1 / (1 / tildeaj + 1 / tildebj + 1 / tildecj + 1 / tildedj)
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
    X2.HBD <- X2.HBD - (sum(a) - sum(tildea))^2 / sum(Var.a)
  }

  STATISTIC <- unname(X2.HBD)
  PARAMETER <- K - 1L

  structure(
    list(
      statistic = c("X-squared" = STATISTIC),
      parameter = c(df = PARAMETER),
      p.value   = pchisq(STATISTIC, PARAMETER, lower.tail = FALSE),
      method    = if (correct)
        "Breslow-Day test for homogeneity of the odds ratios (Tarone corrected)"
      else
        "Breslow-Day test for homogeneity of the odds ratios",
      data.name = DNAME,
      n         = sum(x)
    ),
    class = "htest"
  )
}
