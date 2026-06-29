#' Dunnett's Test for Comparing Several Treatments With a Control
#' 
#' A parametric post hoc test for multiple comparisons of several 
#' treatment groups against a single control group, controlling the 
#' familywise error rate.
#' 
#' Performs Dunnett's test for comparing several treatments with a control.
#' 
#' \code{dunnettTest} does the post hoc pairwise multiple comparisons
#' procedure.
#' 
#' If \code{x} is a list, its elements are taken as the samples to be compared,
#' and hence have to be numeric data vectors.  In this case, \code{g} is
#' ignored, and one can simply use \code{dunnettTest(x)} to perform the test.
#' If the samples are not yet contained in a list, use
#' \code{dunnettTest(list(x, ...))}.
#' 
#' Otherwise, \code{x} must be a numeric data vector, and \code{g} must be a
#' vector or factor object of the same length as \code{x} giving the group for
#' the corresponding elements of \code{x}.
#' 
#' @name dunnettTest
#' @aliases dunnettTest dunnettTest.default dunnettTest.formula
#' @param x a numeric vector of data values, or a list of numeric data vectors.
#' @param g a vector or factor object giving the group for the corresponding
#' elements of \code{x}.  Ignored if \code{x} is a list.
#' @param control the level of the control group against which the others
#' should be tested. Defaults to the first group. Multiple levels are accepted;
#' the calculation is then performed separately for each one. Must be a
#' character string matching one or more levels of \code{g}.
#' @param conf.level confidence level of the interval.
#' @param formula a formula of the form \code{lhs ~ rhs} where \code{lhs} gives
#' the data values and \code{rhs} the corresponding groups.
#' @param data an optional matrix or data frame (or similar: see
#' \code{\link{model.frame}}) containing the variables in the formula
#' \code{formula}.  By default the variables are taken from
#' \code{environment(formula)}.
#' @param subset an optional vector specifying a subset of observations to be
#' used.
#' @param na.action a function which indicates what should happen when the data
#' contain \code{NA}s.  Defaults to \code{getOption("na.action")}.
#' @param \dots further arguments to be passed to or from methods.
#' @return A list of class \code{"PostHocTest"}, containing one matrix named
#' after each control level, with columns \code{diff} giving the difference in
#' observed means, \code{lwr.ci} giving the lower end point of the interval,
#' \code{upr.ci} giving the upper end point, and \code{pval} giving the p-value
#' after adjustment for multiple comparisons.
#' 
#' There are print and plot methods for class \code{"PostHocTest"}. The plot
#' method does not accept \code{xlab}, \code{ylab} or \code{main} arguments and
#' creates its own values for each plot.
#'
#' @seealso
#' \code{\link{print.PostHocTest}},
#' \code{\link{plot.PostHocTest}},
#' \code{\link[mvtnorm]{pmvt}},
#' \code{\link[mvtnorm]{qmvt}}
#' 
#' @references Dunnett C. W. (1955) A multiple comparison procedure for
#' comparing several treatments with a control, \emph{Journal of the American
#' Statistical Association}, 50:1096--1121.
#' 
#' @examples
#' 
#' ## Hollander & Wolfe (1973), 116.
#' ## Mucociliary efficiency from the rate of removal of dust in normal
#' ##  subjects, subjects with obstructive airway disease, and subjects
#' ##  with asbestosis.
#' x <- c(2.9, 3.0, 2.5, 2.6, 3.2) # normal subjects
#' y <- c(3.8, 2.7, 4.0, 2.4)      # with obstructive airway disease
#' z <- c(2.8, 3.4, 3.7, 2.2, 2.0) # with asbestosis
#' 
#' dunnettTest(list(x, y, z))
#' 
#' ## Equivalently,
#' x <- c(x, y, z)
#' g <- factor(rep(1:3, c(5, 4, 5)),
#'             labels = c("Normal subjects",
#'                        "Subjects with obstructive airway disease",
#'                        "Subjects with asbestosis"))
#' 
#' dunnettTest(x, g)
#' 
#' ## Formula interface
#' boxplot(Ozone ~ Month, data = airquality)
#' dunnettTest(Ozone ~ Month, data = airquality)
#' 
#' ## Single control level with adjusted confidence
#' dunnettTest(Ozone ~ Month, data = airquality, control = "8", conf.level = 0.9)
#'
#' ## Multiple control levels
#' dunnettTest(Ozone ~ Month, data = airquality, control = c("5", "8"))
#' 
#' @rdname dunnettTest

#' @family test.posthoc  
#' @concept post-hoc  
#' @concept parametric
#'
#'
#' @export
dunnettTest <- function(x, ...)
  UseMethod("dunnettTest")


# ======================================================================
# Formula method
# ======================================================================

#' @rdname dunnettTest
#' @export
dunnettTest.formula <- function(
    formula,
    data,
    subset,
    na.action,
    ...
) {

  if (missing(formula) || length(formula) != 3L)
    stop("'formula' missing or incorrect")

  ## IMPORTANT!!
  ## --- capture subset / na.action HERE ---
  subset_expr <- if (!missing(subset)) substitute(subset) else NULL
  na_expr     <- if (!missing(na.action)) substitute(na.action) else NULL

  pf <- resolveFormula(
    formula   = formula,
    data      = data,
    subset    = subset_expr,
    na.action = na_expr,
    allowed   = "n-sample-independent"
  )

  rval <- dunnettTest(
    x = pf$x,
    g = pf$group,
    ...
  )

  rval$data.name <- pf$data.name

  rval

}


#' @rdname dunnettTest
#' @export
dunnettTest.default <- function(
    x,
    g,
    control = NULL,
    conf.level = 0.95,
    ...
) {
  
  DG <- resolveGroups(x, g)
  
  x <- DG$x
  g <- DG$groups
  
  N <- DG$n
  k <- DG$k
  
  gn <- DG$group.names
  
  ni <- as.numeric(DG$group.sizes)
  names(ni) <- gn
  
  DNAME <- DG$data.name
  
  if (is.null(control))
    control <- gn[1]
  
  ctrls <- control
  
  if (!all(ctrls %in% gn)) {
    stop(
      gettextf(
        "control level '%s' not found in grouping variable",
        ctrls[!ctrls %in% gn][1]
      )
    )
  }
  
  means <- tapply(x, g, mean)
  
  rss <- tapply(
    x,
    g,
    function(z)
      sum((z - mean(z))^2)
  )
  
  s <- sqrt(sum(rss) / (N - k))
  
  out <- vector("list", length(ctrls))
  
  for (ii in seq_along(ctrls)) {
    
    control <- ctrls[ii]
    
    meandiffs <-
      means[names(means) != control] -
      means[control]
    
    fittedn <- ni[names(ni) != control]
    controln <- ni[control]
    
    Dj <- meandiffs /
      (s * sqrt((1 / fittedn) + (1 / controln)))
    
    Rij <- sqrt(
      fittedn / (fittedn + controln)
    )
    
    R <- outer(Rij, Rij, "*")
    diag(R) <- 1
    
    qvt <- mvtnorm::qmvt(
      p = 1 - (1 - conf.level) / 2,
      df = N - k,
      sigma = R,
      tail = "lower.tail",
      seed = 5L
    )$quantile
    
    lower <-
      meandiffs -
      s * sqrt((1 / fittedn) + (1 / controln)) * qvt
    
    upper <-
      meandiffs +
      s * sqrt((1 / fittedn) + (1 / controln)) * qvt
    
    pval <- numeric(length(meandiffs))
    
    for (i in seq_along(meandiffs)) {
      
      pval[i] <-
        1 -
        mvtnorm::pmvt(
          lower = -abs(Dj[i]),
          upper =  abs(Dj[i]),
          corr  = R,
          delta = rep(0, length(meandiffs)),
          df    = N - k,
          seed = 5L
        )[1]
      
    }
    
    out[[ii]] <- cbind(
      diff   = meandiffs,
      lwr.ci = lower,
      upr.ci = upper,
      pval   = pval
    )
    
    rownames(out[[ii]]) <- paste(
      names(meandiffs),
      control,
      sep = "-"
    )
    
  }
  
  names(out) <- ctrls
  
  class(out) <- "PostHocTest"
  
  attr(out, "conf.level") <- conf.level
  attr(out, "ordered") <- FALSE
  attr(out, "method") <- "Dunnett"
  attr(out, "data.name") <- DNAME
  
  attr(out, "method.str") <- gettextf(
    "\n  Dunnett's test for comparing several treatments with a control : %s\n",
    attr(out, "method")
  )
  
  out
  
}
