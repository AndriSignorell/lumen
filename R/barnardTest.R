
#' Barnard's Unconditional Test for Comparing Two Independent Proportions in a 2 x 2 Table
#'
#' An exact unconditional test for \eqn{2 \times 2}{2x2} contingency
#' tables, offering a more powerful alternative to Fisher's exact test
#' by avoiding conditioning on both marginal totals.
#'
#' Two fundamentally different exact tests exist for the equality of two
#' binomial probabilities: Fisher's exact test (Fisher, 1925), which
#' conditions on both margins, and Barnard's unconditional test (Barnard,
#' 1945), which does not. Fisher's test is by far the more widely used, and
#' Fisher himself objected strongly to Barnard's proposal. For 2 x 2 tables,
#' however, Barnard's test is generally more powerful, as Barnard already
#' pointed out; its limited use is mainly due to its computational cost
#' (Mehta and Senchaudhuri, 2003).
#'
#' **Choice of test.** In most studies only one margin is fixed by design,
#' e.g. the group sizes of a randomised trial, while the number of events is
#' random. Conditioning on both margins, as Fisher's exact test does, then
#' makes the test markedly conservative: because of the discreteness of the
#' conditional distribution, the actual level often lies far below the
#' nominal one, at the cost of power, in particular for the small samples
#' for which exact tests are used. Unconditional tests avoid this. Boschloo's
#' test, which orders the tables by Fisher's p-value, is uniformly more
#' powerful than Fisher's test: its p-value never exceeds Fisher's, while it
#' still keeps the nominal level (Boschloo, 1970). Lydersen et al. (2009)
#' therefore recommend unconditional tests (or Fisher's mid-p) instead of
#' the traditional Fisher test, and the score statistic (`"z-pooled"`) is to
#' be preferred to the Wald statistic (`"z-unpooled"`), which can behave
#' poorly (Mehrotra et al., 2003).
#'
#' For two independent proportions with fixed group sizes, `barnardTest()`
#' is thus the recommended choice; CSM and Boschloo are the most powerful
#' orderings, z-pooled is the most widely used one. [fisher.test()] remains
#' appropriate for r x c tables, for exact conditional inference on the odds
#' ratio, when both margins are fixed by design, and whenever conditional
#' inference is wanted on principle: whether the event margin should be
#' conditioned on is a long-standing debate (Yates, 1984; Little, 1989).
#'
#' **Sampling model.** Unconditional exact tests can be performed for
#' binomial or multinomial models. The binomial model assumes that the row
#' or the column margins (but not both) are fixed by design, the
#' multinomial model only the total sample size. `fixed` selects the model:
#' `1` (default) for fixed rows, `2` for fixed columns, `NA` for the
#' multinomial model. Conditioning on both margins leads to Fisher's exact
#' test, see [fisher.test()].
#'
#' **Binomial model.** The fixed margin defines the two groups, of sizes
#' \eqn{n_1}{n1} and \eqn{n_2}{n2}; the first category of the other
#' classification counts the successes \eqn{x_1}{x1} and \eqn{x_2}{x2}.
#' The tested difference is \eqn{\pi_1 - \pi_2}{pi1 - pi2} with the
#' estimate \eqn{\hat{p}_1-\hat{p}_2}{p1 - p2}, \eqn{\hat{p}_i = x_i/n_i}{pi =
#' xi/ni}, against the null value `delta`. Under the null hypothesis the
#' probability of a table is
#' \deqn{P(x_1, x_2) = \binom{n_1}{x_1}(\pi+\delta)^{x_1}(1-\pi-\delta)^{n_1-x_1}
#'   \binom{n_2}{x_2}\pi^{x_2}(1-\pi)^{n_2-x_2},}{P(x1, x2) =
#'   dbinom(x1, n1, pi + delta) * dbinom(x2, n2, pi),}
#' where \eqn{\pi}{pi} is the unknown nuisance parameter. The tables are
#' ordered by how extreme they are, the p-value for a given \eqn{\pi}{pi}
#' is the probability of all tables at least as extreme as the observed
#' one, and the p-value of the test is its maximum over \eqn{\pi}{pi}. The
#' maximum is searched on a grid of `npNumbers` equally spaced values over
#' the admissible range \eqn{[\max(0, -\delta), \min(1, 1-\delta)]}{[max(0,
#' -delta), min(1, 1 - delta)]}, bounds included, and then refined locally
#' with [optimize()].
#'
#' `method` defines the ordering:
#' \describe{
#'   \item{`"csm"`}{Barnard's original convexity, symmetry and minimum
#'     ordering (CSM; Barnard, 1947). Starting from the most extreme table,
#'     the table that increases the maximal p-value least is added to the
#'     rejection region, step by step. Generally the most powerful choice;
#'     binomial model only.}
#'   \item{`"z-pooled"`}{the score statistic
#'     \deqn{Z = \frac{\hat{p}_1-\hat{p}_2-\delta}{\sqrt{\tilde{p}_1(1-\tilde{p}_1)/n_1 + \tilde{p}_2(1-\tilde{p}_2)/n_2}},}{Z = (p1 - p2 - delta) / sqrt(p1~(1 - p1~)/n1 + p2~(1 - p2~)/n2),}
#'     with the pooled proportion for \eqn{\tilde{p}_1 = \tilde{p}_2}{p1~ =
#'     p2~} if \eqn{\delta = 0}{delta = 0}, and the restricted maximum
#'     likelihood estimates (Miettinen and Nurminen, 1985) otherwise
#'     (Suissa and Shuster, 1985).}
#'   \item{`"z-unpooled"`}{the Wald statistic, with the unrestricted
#'     estimates \eqn{\hat{p}_i}{pi} in the variance (and the same numerator
#'     \eqn{\hat{p}_1-\hat{p}_2-\delta}{p1 - p2 - delta}).}
#'   \item{`"boschloo"`}{the p-value of Fisher's exact test (Boschloo,
#'     1970).}
#'   \item{`"santner-snell"`}{the plain difference \eqn{\hat{p}_1 -
#'     \hat{p}_2}{p1 - p2} (Santner and Snell, 1980).}
#' }
#'
#' Two-sided tests are by default based on the absolute value of the
#' statistic, i.e. on a region that is symmetric about the null value
#' (`tsMethod = "square"`). With `tsMethod = "central"` the two-sided
#' p-value is twice the smaller one-sided p-value.
#'
#' With `npInterval = TRUE` the maximum is only searched within a
#' \eqn{100(1-\beta)}{100(1 - beta)}% Clopper-Pearson interval of the
#' nuisance parameter and \eqn{\beta}{beta} is added to the p-value
#' (Berger and Boos, 1994). The interval is based on the pooled successes,
#' which are binomial only for \eqn{\delta = 0}{delta = 0}; it is
#' therefore not available with a non-zero `delta`, a confidence interval
#' or CSM.
#'
#' The confidence interval for \eqn{\pi_1 - \pi_2}{pi1 - pi2} (binomial
#' model only) is computed if `conf.level` is not `NA`, and is obtained by
#' inverting the test over `delta`. It is computationally demanding, in particular for
#' CSM, where the ordering has to be rebuilt for every candidate value.
#'
#' **Multinomial model.** Only the total \eqn{N} is fixed. The nuisance
#' parameters are the probabilities of the first row and of the first
#' column, and the tested quantity is the difference in the products of
#' the cell proportions, \eqn{(x_{11}x_{22} - x_{12}x_{21})/N^2}{(x11 x22 -
#' x12 x21) / N^2}. Available for the z-pooled, z-unpooled and Boschloo
#' orderings, with `delta = 0` and without confidence interval.
#'
#' If `x` is a matrix, it is taken as a two-dimensional contingency
#' table, and hence its entries should be nonnegative integers. Otherwise,
#' both `x` and `y` must be vectors of the same length. Incomplete cases are
#' removed, the vectors are coerced into factor objects, and the contingency
#' table is computed from these.
#'
#' @param x a numeric vector or a two-dimensional contingency table in
#' matrix form. `x` and `y` can also both be factors.
#' @param y a factor object; ignored if `x` is a matrix.
#' @param alternative a character string specifying the alternative
#' hypothesis, must be one of `"two.sided"` (default), `"greater"` or
#' `"less"`, referring to \eqn{\pi_1 - \pi_2}{pi1 - pi2} relative to
#' `delta`. You can specify just the initial letter.
#' @param method the ordering of the tables, one of `"csm"` (default),
#' `"z-pooled"`, `"z-unpooled"`, `"boschloo"` or `"santner-snell"`. See
#' Details.
#' @param fixed indicates which margin is fixed: `1` for rows (default),
#' `2` for columns, or `NA` for none of both (multinomial model).
#' @param tsMethod the two-sided method, `"square"` (default) or
#' `"central"`. See Details.
#' @param delta the difference \eqn{\pi_1 - \pi_2}{pi1 - pi2} under the null
#' hypothesis, strictly between -1 and 1; non-zero values give tests of
#' non-inferiority or superiority. Binomial model only.
#' @param conf.level confidence level of the interval for
#' \eqn{\pi_1 - \pi_2}{pi1 - pi2} (binomial model only). The default `NA`
#' computes no interval.
#' @param npNumbers the number of grid points for the nuisance parameter,
#' at least 2, default 100.
#' @param npInterval logical, restrict the nuisance parameter to a
#' confidence interval (Berger and Boos, 1994), default `FALSE`.
#' @param beta the error probability of that interval, default 0.001; only
#' used with `npInterval = TRUE`.
#'
#' @return A list with class `"htest"` containing the following
#' components:
#'   \item{`statistic`}{the ordering statistic of the observed table: `z`
#'     (z-pooled, z-unpooled), `d` (Santner-Snell, the difference minus
#'     `delta`) or `p.Fisher` (Boschloo); absent for CSM, whose ordering has
#'     no statistic.}
#'   \item{`parameter`}{the sizes `n1`, `n2` of the two samples (binomial model) or the
#'     total sample size `N` (multinomial model).}
#'   \item{`p.value`}{the p-value of the test.}
#'   \item{`conf.int`}{the confidence interval, if `conf.level` is not
#'     `NA`.}
#'   \item{`estimate`}{the observed difference in proportions (binomial) or
#'     in products of proportions (multinomial).}
#'   \item{`null.value`}{the difference under the null hypothesis.}
#'   \item{`alternative`}{a character string describing the alternative
#'     hypothesis.}
#'   \item{`np`}{the value(s) of the nuisance parameter(s) that maximize the
#'     p-value.}
#'   \item{`np.range`}{the range of nuisance parameters considered.}
#'   \item{`model`, `tsmethod`}{the sampling model and the two-sided method.}
#'   \item{`method`}{a character string describing the test.}
#'   \item{`data.name`}{a character string giving the name of the data.}
#'
#' @note
#' The algorithms follow the \pkg{Exact} package by Peter Calhoun
#' (`Exact::exact.test()`), with these deliberate differences:
#' \itemize{
#'   \item the rejection region of the ordering statistics consists of
#'     exactly the tables that are at least as extreme as the observed one.
#'     Exact builds it by walking the boundary and filling rows and columns
#'     by convexity, which in two-sided tests can admit tables that are
#'     less extreme (e.g. `matrix(c(3, 2, 1, 0), 2)`, z-pooled: 0.783
#'     instead of 0.568);
#'   \item the nuisance grid includes the bounds of the admissible range
#'     (Exact uses 0.00001 and 0.99999), where the maximum can lie;
#'   \item the z-unpooled statistic subtracts `delta` in the numerator, as
#'     the z-pooled one does;
#'   \item the stored CSM orderings of the companion package \pkg{ExactData}
#'     are not used: the compiled ordering is fast enough to be computed on
#'     the fly, for every `npNumbers`;
#'   \item when the grid maximum is attained at several nuisance values, the
#'     refined p-value is the largest of the refined maxima (Exact takes the
#'     smallest).
#' }
#'
#' @seealso [fisher.test()], [prop.test()]
#'
#' @references Barnard, G.A. (1945) A new test for 2x2 tables.
#' *Nature*, 156:177.
#'
#' Barnard, G.A. (1947) Significance tests for 2x2 tables.
#' *Biometrika*, 34:123-138.
#'
#' Berger, R.L. and Boos, D.D. (1994) P values maximised over a confidence
#' set for the nuisance parameter. *Journal of the American Statistical
#' Association*, 89:1012-1016.
#'
#' Boschloo, R.D. (1970) Raised conditional level of significance for the
#' 2x2-table when testing the equality of two probabilities. *Statistica
#' Neerlandica*, 24:1-35.
#'
#' Calhoun, P. (2024) Exact: Unconditional Exact Test. R package version
#' 3.3.
#'
#' Little, R.J.A. (1989) Testing the equality of two independent binomial
#' proportions. *The American Statistician*, 43:283-288.
#'
#' Lydersen, S., Fagerland, M.W. and Laake, P. (2009) Recommended tests for
#' association in 2x2 tables. *Statistics in Medicine*, 28:1159-1175.
#'
#' Mehrotra, D.V., Chan, I.S.F. and Berger, R.L. (2003) A cautionary note on
#' exact unconditional inference for a difference between two independent
#' binomial proportions. *Biometrics*, 59:441-450.
#'
#' Mehta, C.R., Senchaudhuri, P. (2003) Conditional versus unconditional
#' exact tests for comparing two binomials. Cytel Software Corporation,
#' technical report.
#'
#' Miettinen, O. and Nurminen, M. (1985) Comparative analysis of two rates.
#' *Statistics in Medicine*, 4:213-226.
#'
#' Santner, T.J. and Snell, M.K. (1980) Small-sample confidence intervals
#' for p1 - p2 and p1/p2 in 2x2 contingency tables. *Journal of the
#' American Statistical Association*, 75:386-394.
#'
#' Suissa, S. and Shuster, J. J. (1985) Exact unconditional sample sizes for
#' the 2x2 binomial trial, *Journal of the Royal Statistical Society*,
#' Ser. A, 148, 317-327.
#'
#' Yates, F. (1984) Tests of significance for 2 x 2 contingency tables.
#' *Journal of the Royal Statistical Society*, Ser. A, 147:426-463.
#'
#' @examples
#' tab <- as.table(matrix(c(8, 14, 1, 3), nrow=2,
#'                 dimnames=list(treat=c("I","II"), out=c("I","II"))))
#' barnardTest(tab, method="z-pooled")
#'
#' # the default CSM ordering
#' barnardTest(tab)
#'
#' # Mehta and Senchaudhuri (2003), vaccine example
#' tab <- as.table(matrix(c(7, 12, 8, 3), nrow=2,
#'                        dimnames=list(treat=c("vaccine","placebo"),
#'                                      infection=c("yes","no"))))
#' barnardTest(tab, alternative="less")
#'
#' # non-inferiority margin of -0.1, with confidence interval
#' barnardTest(tab, method="z-pooled", alternative="greater", delta=-0.1,
#'             conf.level=0.95)
#'
#' # multinomial model: only the total is fixed
#' barnardTest(tab, method="boschloo", fixed=NA)
#'
#' @family test.categorical
#' @concept categorical-test
#' @concept exact-test
#'
#' @export
barnardTest <- function(x,
                        y = NULL,
                        alternative = c("two.sided", "less", "greater"),
                        method = c("csm", "z-pooled", "z-unpooled",
                                   "boschloo", "santner-snell"),
                        fixed = 1,
                        tsMethod = c("square", "central"),
                        delta = 0,
                        conf.level = NA,
                        npNumbers = 100,
                        npInterval = FALSE,
                        beta = 0.001) {

  alternative <- match.arg(alternative)
  method      <- match.arg(method)
  tsMethod    <- match.arg(tsMethod)

  # deparse1: deparse() splits long expressions into several strings
  DNAME <- if (is.null(y) || !is.null(dim(x))) deparse1(substitute(x)) else
    paste(deparse1(substitute(x)), "and", deparse1(substitute(y)))

  # --- validate the arguments before touching the data -------------------
  if (length(fixed) == 2 && setequal(fixed, 1:2))
    stop("use fisher.test() if both margins, rows AND columns, are fixed")
  if (length(fixed) != 1 ||
      (!is.na(fixed) && (!is.numeric(fixed) || !fixed %in% c(1, 2))))
    stop("'fixed' must be 1 (rows), 2 (columns) or NA (none of both)")

  withCI <- !(length(conf.level) == 1L && is.na(conf.level))
  if (withCI) checkConfLevel(conf.level)
  checkFlag(npInterval)
  if (!is.numeric(delta) || length(delta) != 1L || !is.finite(delta) ||
      delta <= -1 || delta >= 1)
    stop("'delta' must be a single number strictly between -1 and 1")
  if (!is.numeric(npNumbers) || length(npNumbers) != 1L ||
      !is.finite(npNumbers) || npNumbers < 2 || npNumbers != round(npNumbers))
    stop("'npNumbers' must be an integer of at least 2")
  npNumbers <- as.integer(npNumbers)
  if (npInterval && (!is.numeric(beta) || length(beta) != 1L ||
                     !is.finite(beta) || beta < 0 || beta > 1))
    stop("'beta' must be a single number in [0, 1]")

  multinomial <- is.na(fixed)

  if (multinomial) {
    if (method %in% c("csm", "santner-snell"))
      stop(gettextf("the %s ordering is not available for the multinomial model (fixed = NA); choose another 'method'",
                    sQuote(method)), call. = FALSE)
    if (delta != 0 || withCI)
      stop("non-zero 'delta' and confidence intervals are only available ",
           "for the binomial model", call. = FALSE)
  }
  if (npInterval && method == "csm")
    stop("'npInterval' cannot be used with the CSM ordering", call. = FALSE)
  # the pooled successes are binomial only under pi1 = pi2
  if (npInterval && (delta != 0 || withCI))
    stop("'npInterval' requires delta = 0 and no confidence interval",
         call. = FALSE)
  # Fisher's p-value does not depend on delta, so a two-sided region based
  # on it cannot be centred at a non-zero null value
  if (alternative == "two.sided" && tsMethod == "square" &&
      method == "boschloo" && (delta != 0 || withCI))
    stop("a two-sided Boschloo test with non-zero 'delta' or a confidence interval needs tsMethod = \"central\"",
         call. = FALSE)

  # --- the table ----------------------------------------------------------
  CT <- resolveContingency(x, y, integerCounts = FALSE, dataName = DNAME)
  tab <- CT$table

  if (nrow(tab) != 2L || ncol(tab) != 2L)
    stop("'x' must be a 2 x 2 contingency table")
  if (any(tab != round(tab)))
    stop("all entries of 'x' must be integer counts")
  tab <- matrix(as.integer(tab), 2L, 2L)

  methodName <- switch(method,
                       "csm"           = "CSM ordering",
                       "z-pooled"      = "z-pooled ordering",
                       "z-unpooled"    = "z-unpooled ordering",
                       "boschloo"      = "Boschloo ordering",
                       "santner-snell" = "Santner-Snell ordering")
  # name of the ordering statistic (none for CSM)
  STATNAME <- switch(method,
                     "z-pooled" = , "z-unpooled" = "z",
                     "santner-snell" = "d",
                     "boschloo" = "p.Fisher",
                     "csm" = NA_character_)

  METHOD <- paste0("Barnard's unconditional exact test (", methodName,
                   if (npInterval) ", interval approach", ")")

  if (multinomial) {

    res <- .barnardMultinomial(tab, alternative = alternative,
                               method = method, npNumbers = npNumbers,
                               npInterval = npInterval, beta = beta)
    N <- sum(tab)
    ESTIMATE  <- (tab[1, 1] * tab[2, 2] - tab[1, 2] * tab[2, 1]) / N^2
    PARAMETER <- c(N = N)
    NULLVALUE <- 0
    names(ESTIMATE) <- names(NULLVALUE) <- "difference in product of proportion"
    NP <- list("first nuisance parameter"  = res$np1,
               "second nuisance parameter" = res$np2)
    NPRANGE <- list("first nuisance parameter range"  = res$np1.range,
                    "second nuisance parameter range" = res$np2.range)
    CINT <- NULL

  } else {

    # groups as columns, successes in the first row
    if (fixed == 1) tab <- t(tab)
    Ns <- colSums(tab)
    if (any(Ns == 0L))
      stop("both groups must contain at least one observation", call. = FALSE)

    res <- .barnardBinomial(tab, alternative = alternative, method = method,
                            tsMethod = tsMethod, delta = delta,
                            npNumbers = npNumbers, npInterval = npInterval,
                            beta = beta)

    ESTIMATE  <- tab[1, 1] / Ns[1] - tab[1, 2] / Ns[2]
    PARAMETER <- setNames(Ns, c("n1", "n2"))
    NULLVALUE <- delta
    names(ESTIMATE) <- names(NULLVALUE) <- "difference in proportion"
    NP <- list("nuisance parameter" = res$np)
    NPRANGE <- list("nuisance parameter range" = res$np.range)

    CINT <- if (withCI) {
      ci <- .barnardConfInt(tab, conf.level = conf.level,
                            estimate = unname(ESTIMATE),
                            alternative = alternative, method = method,
                            tsMethod = tsMethod, npNumbers = npNumbers,
                            npInterval = npInterval, beta = beta)
      structure(ci, conf.level = conf.level)
    }
  }

  structure(
    # CSM orders the tables without a statistic: no 'statistic' component,
    # so that print.htest() does not show "test statistic = NA"
    list(statistic   = if (!is.na(res$statistic))
                         setNames(res$statistic, STATNAME),
         parameter   = PARAMETER,
         p.value     = res$p.value,
         conf.int    = CINT,
         estimate    = ESTIMATE,
         null.value  = NULLVALUE,
         alternative = alternative,
         np          = NP,
         np.range    = NPRANGE,
         model       = if (multinomial) "multinomial" else "binomial",
         method      = METHOD,
         tsmethod    = tsMethod,
         data.name   = DNAME),
    class = "htest")
}



# == internal helpers ===================================================

# Clopper-Pearson interval, as Exact's binom.CI()
.barnardClopperPearson <- function(x, n, conf.level) {
  alpha <- (1 - conf.level) / 2
  c(if (x == 0) 0 else qbeta(alpha, x, n - x + 1),
    if (x == n) 1 else qbeta(1 - alpha, x + 1, n - x))
}


# Nuisance grid: the admissible range of pi (the success probability of the
# second group, pi + delta that of the first) including its bounds, where
# the maximum may lie; with npInterval (delta = 0 only) restricted to a
# Clopper-Pearson interval of the pooled success probability
.barnardGrid <- function(successes, N, delta, npNumbers, npInterval, beta) {
  lo <- max(0, -delta)
  hi <- min(1, 1 - delta)
  if (npInterval && beta != 0) {
    ci <- .barnardClopperPearson(successes, N, conf.level = 1 - beta)
    lo <- max(lo, ci[1L])
    hi <- min(hi, ci[2L])
  }
  seq(lo, hi, length.out = npNumbers)
}


# p-value of the binomial model; 'tab' has the groups as columns and the
# successes in the first row
.barnardBinomial <- function(tab, alternative, method, tsMethod, delta,
                             npNumbers, npInterval, beta) {

  n1 <- sum(tab[, 1L]); n2 <- sum(tab[, 2L])
  x  <- tab[1L, 1L];    y  <- tab[1L, 2L]

  # orient the table so that the observed difference lies on the side of
  # the "less" alternative
  swap <- alternative == "greater" ||
    (alternative == "two.sided" && (x / n1 - y / n2) > delta)
  if (swap) {
    tmp <- x;  x  <- y;  y  <- tmp
    tmp <- n1; n1 <- n2; n2 <- tmp
    delta <- -delta
    if (alternative == "greater") alternative <- "less"
  }

  doubleP <- alternative == "two.sided" && tsMethod == "central"
  if (doubleP) alternative <- "less"
  twoSided <- alternative == "two.sided"

  grid <- .barnardGrid(x + y, n1 + n2, delta, npNumbers, npInterval, beta)
  if (!npInterval) beta <- 0

  if (method == "csm") {
    region <- .barnardCSM(n1, n2, x, y, grid, delta, twoSided, doubleP)
    TXO <- NA_real_
  } else {
    # every table is compared with the observed one. (Exact walks the
    # boundary row by row and fills rows and columns by convexity, which
    # admits tables that are not as extreme in two-sided tests.)
    TX <- .barnardStat(n1, n2, method, alternative, delta)
    TXO <- TX[x + 1L, y + 1L]
    if (anyNA(TX))
      stop("test statistic is NA", call. = FALSE)
    region <- if (twoSided && method != "boschloo") abs(TX) >= abs(TXO)
              else TX <= TXO
    storage.mode(region) <- "integer"
  }

  prob <- signif((1 + doubleP) * .barnardProb(region, grid, delta) + beta, 12)
  prob[!is.na(prob) & prob > 1] <- 1
  pvalue <- max(prob, na.rm = TRUE)
  np <- grid[which(prob == pvalue)]

  # refine the maximum locally, over the cells next to each grid maximum
  # (a local search: it improves on the grid but cannot guarantee the
  # global supremum)
  step <- diff(grid[1:2])
  ref <- lapply(np, function(p0)
    optimize(function(p) .barnardProb(region, p, delta),
             interval = c(max(grid[1L], p0 - step),
                          min(grid[npNumbers], p0 + step)),
             maximum = TRUE))
  refP <- signif((1 + doubleP) *
                   vapply(ref, `[[`, numeric(1), "objective") + beta, 12)
  if (any(!is.na(refP)) && max(refP, na.rm = TRUE) > pvalue) {
    pvalue <- min(1, max(refP, na.rm = TRUE))
    # symmetric problems have tied maxima (p and 1 - p): keep all of them
    np <- unique(vapply(ref, `[[`, numeric(1), "maximum")[
      !is.na(refP) & refP == max(refP, na.rm = TRUE)])
  }

  if (swap && method %in% c("z-pooled", "z-unpooled", "santner-snell"))
    TXO <- -TXO

  # the nuisance parameter is the success probability of the second group;
  # after the swap that is the internal first group, pi + delta
  shift <- if (swap) delta else 0

  list(p.value = pvalue, statistic = TXO, np = np + shift,
       np.range = range(grid) + shift)
}


# p-value of the multinomial model (Exact's multinomCode() and
# maxPvalueMultinom())
.barnardMultinomial <- function(tab, alternative, method, npNumbers,
                                npInterval, beta) {

  N <- sum(tab)
  mm <- .barnardMultinom(N, tab[1L, 1L], tab[1L, 2L], tab[2L, 1L],
                         method, alternative)

  if (npInterval) {
    int1 <- .barnardGrid(sum(tab[1L, ]), N, 0, npNumbers, TRUE, beta)
    int2 <- .barnardGrid(sum(tab[, 1L]), N, 0, npNumbers, TRUE, beta)
  } else {
    int1 <- int2 <- seq(0, 1, length.out = npNumbers)
    beta <- 0
  }

  # (p1, p2) and the swapped pairs, restricted to p1 <= p2
  pairs <- expand.grid(p1 = int1, p2 = int2)
  pairs <- rbind(pairs, setNames(rev(pairs), names(pairs)))
  pairs <- pairs[pairs$p1 <= pairs$p2 & !duplicated(pairs), ]

  u1 <- unique(pairs$p1)
  u2 <- unique(pairs$p2)
  P  <- .barnardMultinomProb(mm$W, u1, u2)

  maxProb <- np1 <- numeric(length(u2))
  for (k in seq_along(u2)) {
    p1k <- pairs$p1[pairs$p2 == u2[k]]
    vals <- P[match(p1k, u1), k]
    maxProb[k] <- max(vals, na.rm = TRUE)
    np1[k] <- p1k[vals == maxProb[k]][1L]
  }
  maxProb <- signif(maxProb + beta, 12)
  maxProb[!is.na(maxProb) & maxProb > 1] <- 1

  best <- maxProb == max(maxProb)
  reverse <- np1[best] %in% int1
  pvalue <- max(maxProb)
  NP1 <- c(np1[best], u2[best])[c(reverse, !reverse)]
  NP2 <- c(np1[best], u2[best])[c(!reverse, reverse)]

  # refine the grid maxima locally, as in the binomial model (Exact uses
  # the grid value, which can fall short of the supremum). The best cell
  # centre is a further start: on a coarse grid all grid points can lie
  # where the probability and its gradient vanish (e.g. diag(2) with
  # npNumbers = 2: zero at the corners, 0.25 at the centre). A local search
  # improves on the grid, it does not guarantee the global supremum.
  step1 <- diff(int1[1:2])
  step2 <- diff(int2[1:2])
  c1 <- (int1[-1L] + int1[-npNumbers]) / 2
  c2 <- (int2[-1L] + int2[-npNumbers]) / 2
  Pc <- .barnardMultinomProb(mm$W, c1, c2)
  bc <- which(Pc == max(Pc), arr.ind = TRUE)[1L, ]
  NP1s <- c(NP1, c1[bc[1L]])
  NP2s <- c(NP2, c2[bc[2L]])

  ref <- Map(function(a, b) {
    lo <- c(max(min(int1), a - step1), max(min(int2), b - step2))
    hi <- c(min(max(int1), a + step1), min(max(int2), b + step2))
    o <- optim(c(a, b), function(p)
      -.barnardMultinomProb(mm$W, p[1L], p[2L])[1L, 1L],
      method = "L-BFGS-B", lower = lo, upper = hi)
    c(o$par, signif(-o$value + beta, 12))
  }, NP1s, NP2s)
  ref <- do.call(rbind, ref)
  if (max(ref[, 3L]) > pvalue) {
    keep <- ref[, 3L] == max(ref[, 3L])
    pvalue <- min(1, max(ref[, 3L]))
    NP1 <- ref[keep, 1L]
    NP2 <- ref[keep, 2L]
  }

  list(p.value   = pvalue,
       statistic = mm$statistic,
       np1 = NP1, np2 = NP2,
       np1.range = range(int1), np2.range = range(int2))
}


# all roots of f on a grid of n subintervals (as rootSolve::uniroot.all)
.unirootAll <- function(f, lower, upper, n = 100, tol = 1e-4) {
  xseq <- seq(lower, upper, length.out = n + 1L)
  mod  <- f(xseq)
  roots <- xseq[which(mod == 0)]
  for (i in which(mod[-(n + 1L)] * mod[-1L] < 0))
    roots <- c(roots, uniroot(f, lower = xseq[i], upper = xseq[i + 1L],
                              tol = tol)$root)
  roots
}


# confidence interval by test inversion over delta (Exact's confInt())
.barnardConfInt <- function(tab, conf.level, estimate, alternative, method,
                            tsMethod, npNumbers, npInterval, beta) {

  alpha <- 1 - conf.level
  rootFun <- function(d)
    vapply(d, function(dd)
      alpha - .barnardBinomial(tab, alternative = alternative,
                               method = method, tsMethod = tsMethod,
                               delta = dd, npNumbers = npNumbers,
                               npInterval = npInterval,
                               beta = beta)$p.value, numeric(1))

  # after the largest (smallest) root, search beyond it once more, as a
  # p-value function need not be monotone in delta
  UB <- 1
  if (alternative != "greater" && estimate < 0.9998) {
    prev <- 999; lower <- max(estimate, -0.9998)
    while (UB != prev && lower < 0.9998) {
      prev <- UB
      r <- .unirootAll(rootFun, lower, 0.9998)
      if (length(r)) UB <- max(r)
      lower <- UB + 0.001
    }
    if (UB >= 0.9998) UB <- 1
  }

  LB <- -1
  if (alternative != "less" && estimate > -0.9998) {
    prev <- 999; upper <- min(estimate, 0.9998)
    while (LB != prev && upper > -0.9998) {
      prev <- LB
      r <- .unirootAll(rootFun, -0.9998, upper)
      if (length(r)) LB <- min(r)
      upper <- LB - 0.001
    }
    if (LB <= -0.9998) LB <- -1
  }

  c(LB, UB)
}
