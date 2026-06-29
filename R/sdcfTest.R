
#' Dwass-Steel-Critchlow-Fligner All-Pairs Test
#'
#' Performs the Dwass-Steel-Critchlow-Fligner (DSCF) all-pairs
#' multiple-comparison procedure for independent samples.
#'
#' The DSCF test is a nonparametric all-pairs comparison procedure based
#' on pairwise Wilcoxon rank-sum statistics and provides strong control
#' of the family-wise error rate. It may be viewed as a nonparametric
#' analogue of Tukey's all-pairs procedure and is generally more
#' powerful than the classical Nemenyi test.
#'
#' For each pair of groups, observations are ranked jointly within the
#' two groups being compared, yielding a pairwise Mann-Whitney statistic.
#' Test statistics are transformed according to the
#' Dwass-Steel-Critchlow-Fligner procedure and evaluated using the
#' studentized range distribution.
#'
#' If \code{x} is a list, its elements are taken as the samples to be
#' compared and hence have to be numeric data vectors. In this case,
#' \code{g} is ignored and one can simply use \code{sdcfTest(x)}.
#'
#' Otherwise, \code{x} must be a numeric vector and \code{g} a grouping
#' factor (or vector coercible to a factor) of the same length.
#'
#' @name sdcfTest
#' @aliases sdcfTest sdcfTest.default sdcfTest.formula
#'
#' @param x a numeric vector of observations, or a list of numeric
#'   vectors.
#' @param g a grouping factor corresponding to \code{x}. Ignored if
#'   \code{x} is a list.
#' @param output character string specifying the output format. One of
#'   \code{"list"} (default) or \code{"matrix"}.
#' @param formula a formula of the form \code{response ~ group}.
#' @param data an optional data frame containing the variables in
#'   \code{formula}.
#' @param subset an optional expression specifying a subset of
#'   observations to be used.
#' @param na.action a function specifying how missing values should be
#'   handled.
#' @param \dots further arguments passed to methods.
#'
#' @return An object of class \code{c("rankTest", "htest")} containing:
#' \describe{
#'   \item{res}{
#'     comparison results. For \code{output="list"} a matrix with
#'     columns \code{z} and \code{pval}; for
#'     \code{output="matrix"} a symmetric matrix of adjusted p-values
#'     with diagonal 1.
#'   }
#'   \item{pmat}{
#'     symmetric matrix of adjusted p-values with diagonal 1.
#'   }
#' }
#'
#' Additional information is stored in attributes:
#' \code{method}, \code{output}, \code{main}, and
#' \code{data.name}.
#'
#' @details
#' The DSCF procedure performs all pairwise comparisons using
#' pair-specific rankings. For each pair of groups, a Mann-Whitney-type
#' statistic is computed and standardized using a tie-corrected variance
#' estimate. Adjusted p-values are obtained from the studentized range
#' distribution and control the family-wise error rate across all
#' pairwise comparisons.
#'
#' The implementation reproduces the procedure described by
#' Dwass (1960), Steel (1960), Critchlow and Fligner (1991), and is
#' equivalent to the implementation in
#' \code{PMCMRplus::dscfAllPairsTest()}.
#'
#' @references
#'
#' Dwass, M. (1960).
#' Some k-sample rank-order tests.
#' In: \emph{Contributions to Probability and Statistics},
#' Stanford University Press, 198--202.
#'
#' Steel, R. G. D. (1960).
#' A rank sum test for comparing all pairs of treatments.
#' \emph{Technometrics}, \strong{2}, 197--207.
#'
#' Critchlow, D. E., & Fligner, M. A. (1991).
#' On distribution-free multiple comparisons in the one-way analysis of
#' variance.
#' \emph{Communications in Statistics - Theory and Methods},
#' \strong{20}, 127--139.
#'
#' @seealso
#' \code{\link{kruskal.test}}
#'
#' @examples
#' ## Hollander & Wolfe style example
#' x <- c(2.9, 3.0, 2.5, 2.6, 3.2)
#' y <- c(3.8, 2.7, 4.0, 2.4)
#' z <- c(2.8, 3.4, 3.7, 2.2, 2.0)
#'
#' sdcfTest(list(x, y, z))
#'
#' x <- c(x, y, z)
#' g <- factor(rep(1:3, c(5, 4, 5)))
#'
#' sdcfTest(x, g)
#'
#' ## Matrix output
#' sdcfTest(x, g, output = "matrix")
#'
#' ## Formula interface
#' sdcfTest(Ozone ~ Month, data = airquality)
#'
#' @family test.posthoc
#' @concept multiple-comparisons
#' @concept nonparametric
#' @concept hypothesis-testing


#' @export
sdcfTest <- function(x, ...)
  UseMethod("sdcfTest")


# ======================================================================
# Formula method
# ======================================================================

#' @rdname sdcfTest
#' @export
sdcfTest.formula <- function(
    formula,
    data,
    subset,
    na.action,
    ...
) {
  
  if (missing(formula) ||
      (length(formula) != 3L) ||
      (length(attr(terms(formula[-2L]), "term.labels")) != 1L))
    stop("'formula' missing or incorrect")
  
  subset_expr <- if (!missing(subset))
    substitute(subset)
  else
    NULL
  
  na_expr <- if (!missing(na.action))
    substitute(na.action)
  else
    NULL
  
  pf <- resolveFormula(
    formula   = formula,
    data      = data,
    subset    = subset_expr,
    na.action = na_expr,
    allowed   = "n-sample-independent"
  )
  
  out <- sdcfTest(
    x = pf$x,
    g = pf$group,
    ...
  )
  
  # consistent with dunnTest / conoverTest pattern
  out$data.name <- pf$data.name
  
  out
}


#' @rdname sdcfTest
#' @export
sdcfTest.default <- function(
    x,
    g,
    output = c(
      "list",
      "matrix"
    ),
    ...
) {
  
  output <- match.arg(output)
  
  DG <- resolveGroups(x, g)
  
  x  <- DG$x
  g  <- DG$groups
  
  gn <- DG$group.names
  k  <- DG$k
  
  comb  <- utils::combn(seq_len(k), 2L)
  npair <- ncol(comb)
  
  z    <- numeric(npair)
  pval <- numeric(npair)
  lbl  <- character(npair)
  
  for (ii in seq_len(npair)) {
    
    i <- comb[1L, ii]
    j <- comb[2L, ii]
    
    xi <- x[g == gn[i]]
    xj <- x[g == gn[j]]
    
    ni <- length(xi)
    nj <- length(xj)
    
    xx <- c(xi, xj)
    
    r <- rank(xx, ties.method = "average")
    
    Ri <- sum(r[seq_len(ni)])
    Rj <- sum(r[ni + seq_len(nj)])
    
    U <- c(
      ni * nj + ni * (ni + 1L) / 2 - Ri,
      ni * nj + nj * (nj + 1L) / 2 - Rj
    )
    
    Umn <- min(U)
    S   <- ni + nj
    
    tt <- table(r)
    C  <- sum((tt^3 - tt) / 12)
    
    VAR <-
      (ni * nj / (S * (S - 1L))) *
      ((S^3 - S) / 12 - C)
    
    if (VAR <= 0) {
      z[ii]    <- 0
      pval[ii] <- 1
    } else {
      z[ii] <-
        sqrt(2) *
        (Umn - ni * nj / 2) /
        sqrt(VAR)
      
      pval[ii] <-
        ptukey(
          abs(z[ii]),
          nmeans = k,
          df     = Inf,
          lower.tail = FALSE
        )
    }
    
    lbl[ii] <- paste(gn[i], gn[j], sep = "-")
  }
  
  # --- result table -------------------------------------------------------
  
  res.mat <- cbind(z = z, pval = pval)
  rownames(res.mat) <- lbl
  
  # --- p-value matrix -----------------------------------------------------
  
  pmat <- matrix(
    NA_real_,
    nrow = k,
    ncol = k,
    dimnames = list(gn, gn)
  )
  
  pmat[lower.tri(pmat)] <- pval
  
  pmatxt <- pmat
  pmatxt[upper.tri(pmatxt)] <- t(pmatxt)[upper.tri(pmatxt)]
  diag(pmatxt) <- 1
  
  attr(pmatxt, "lbl") <- apply(
    pmatxt,
    1,
    function(z)
      paste(rownames(pmatxt)[!is.na(z) & z < 0.05], collapse = ",")
  )
  
  # --- output -------------------------------------------------------------
  
  out <- list()
  
  if (output == "list") {
    out$res <- res.mat
  } else {
    # matrix output: return the full symmetric p-value matrix (diagonal = 1)
    out$res <- pmatxt
  }
  
  out$pmat <- pmatxt
  
  class(out) <- c("rankTest", "htest")
  
  attr(out, "main")      <- "Steel-Dwass-Critchlow-Fligner all-pairs test"
  attr(out, "method")    <- "Dwass-Steel-Critchlow-Fligner"
  attr(out, "output")    <- output
  attr(out, "data.name") <- DG$data.name
  
  out
}

