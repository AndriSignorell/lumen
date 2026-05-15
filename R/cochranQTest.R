
#' Cochran's Q Test
#' 
#' A nonparametric test for dependent samples with dichotomous data, 
#' assessing whether proportions differ across multiple conditions or 
#' time points. It generalizes the McNemar test to more than two groups.
#' 
#' Performs Cochran's Q test for related samples with a binary response.
#' The test is appropriate for unreplicated complete block designs (i.e.,
#' matched or paired data), where each block contains exactly one observation
#' for each group.
#'
#' Cochran's Q test is a nonparametric method for testing whether the
#' proportions of a binary outcome differ across multiple related groups,
#' while accounting for block effects. It can be regarded as an extension
#' of McNemar's test to more than two groups.
#'
#' The null hypothesis is that, apart from block effects, the probability
#' of a "success" (coded as 1) is the same in all groups.
#'
#' If \code{y} is a matrix, groups and blocks are inferred from the column
#' and row indices, respectively. Missing values in \code{y} lead to the
#' removal of the corresponding blocks. Missing values are not allowed in
#' \code{groups} or \code{blocks}.
#'
#' The asymptotic test statistic is computed directly. Optionally, exact or
#' Monte Carlo inference can be obtained via the \pkg{coin} package.
#'
#' Cochran's Q test is closely related to the Friedman test, but is specifically
#' designed for binary (0/1) responses. 
#' 
#' 
#' @name cochranQTest
#' @aliases cochranQTest cochranQTest.default cochranQTest.formula
#'
#' @param y either a numeric vector of data values, or a data matrix.
#' @param method "asymptotic" (default), or "approximate"
#' @param nresample number of Monte Carlo replicates (for approximate)
#' @param groups a vector giving the group for the corresponding elements of y
#' if this is a vector; ignored if y is a matrix. If not a factor object, it is
#' coerced to one.
#' @param blocks a vector giving the block for the corresponding elements of y
#' if this is a vector; ignored if y is a matrix. If not a factor object, it is
#' coerced to one.
#' @param formula a formula of the form \code{y ~ groups | blocks}.
#' @param data an optional matrix or data frame (or similar: see
#' \code{\link{model.frame}}) containing the variables in the formula. By
#' default the variables are taken from \code{environment(formula)}.
#' @param subset an optional vector specifying a subset of observations to be
#' used.
#' @param na.action a function which indicates what should happen when the data
#' contain \code{NA}s. Defaults to \code{getOption("na.action")}.
#' @param \dots further arguments to be passed to or from methods.
#' @return
#' 
#' A list with class \code{htest} containing the following components:
#' 
#' \item{statistic}{the value of Cochran's chi-squared statistic.}
#' \item{parameter}{the degrees of freedom of the approximate chi-squared
#' distribution of the test statistic.} \item{p.value}{the p-value of the
#' test.} \item{method}{the character string "Cochran's Q-Test".}
#' \item{data.name}{a character string giving the names of the data.}
#' 
#' @references Cochran, W.G. (1950) The Comparison of Percentages in Matched
#' Samples. \emph{Biometrika}. 37 (3/4): 256-266.
#' doi:10.1093/biomet/37.3-4.256. JSTOR 2332378.
#' 
#' @examples
#' 
#' # example in: 
#' # http://support.sas.com/documentation/cdl/en/statugfreq/63124/PDF/default/statugfreq.pdf
#' # pp. S. 1824
#' 
#' # create the dataset
#' d.frm <- expand.grid(A=c("F","U"), B=c("F","U"), C=c("F","U"))[
#'             rep(1:8, c(6,2,2,6,16,4,4,6)), ]
#' row.names(d.frm) <- NULL
#' 
#' # rearrange to long shape    
#' d.long <- reshape(d.frm, varying=1:3, times=names(d.frm)[c(1:3)], 
#'                   v.names="resp", direction="long")
#' 
#' 
#' # after having done the hard work of data organisation, performing the test is a piece of cake....
#' cochranQTest(resp ~ time | id, data=d.long)
#' 
#' # and let's perform a post hoc analysis using mcnemar's test
#' z <- split(d.long, f=d.long$time)
#' pairwise.table(function(i, j) { 
#'     mcnemar.test(z[[i]]$resp, z[[j]]$resp, correct=FALSE)$p.value
#'   }, 
#'   level.names = names(z), 
#'   p.adjust.method = "fdr"
#' )
#' 




#' @rdname cochranQTest
#' @family test.goodnessfit
#' @concept hypothesis-testing
#' @concept nonparametric
#'
#'
#' @export
cochranQTest <- function(y, ...) {
  UseMethod("cochranQTest")
}



# =========================
# DEFAULT METHOD
# =========================

#' @rdname cochranQTest
#' @export
cochranQTest.default <- function(y, groups, blocks,
                                 method    = c("asymptotic", "approximate"),
                                 nresample = 1e4,
                                 na.action = na.omit,
                                 ...) {
  
  method <- match.arg(method)
  DNAME  <- deparse1(substitute(y))
  
  if (is.matrix(y)) {
    
    # validate binary
    if (!all(y %in% c(0L, 1L, NA)))
      stop("'y' must contain only 0/1 values")
    
    groups <- factor(c(col(y)))
    blocks <- factor(c(row(y)))
    
  } else {
    
    if (anyNA(groups) || anyNA(blocks))
      stop("NA's are not allowed in 'groups' or 'blocks'")
    
    if (any(diff(c(length(y), length(groups), length(blocks))) != 0L))
      stop("'y', 'groups' and 'blocks' must have the same length")
    
    DNAME <- paste0(DNAME, ", ", deparse1(substitute(groups)),
                    " and ", deparse1(substitute(blocks)))
    
    if (any(table(groups, blocks) != 1L))
      stop("not an unreplicated complete block design")
    
    groups <- factor(groups)
    blocks <- factor(blocks)
    o      <- order(groups, blocks)
    y      <- asBinary(y[o])
    groups <- groups[o]
    blocks <- blocks[o]
  }
  
  ## -------------------------------------------------------------------
  ## Approximate via coin
  ## -------------------------------------------------------------------
  if (method == "approximate") {
    
    if (!requireNamespace("coin", quietly = TRUE))
      stop("package 'coin' required for method = 'approximate'")
    
    df <- data.frame(y = factor(y), groups = groups, blocks = blocks)
    
    res <- coin::independence_test(
      y ~ groups | blocks,
      data         = df,
      teststat     = "quad",
      distribution = coin::approximate(nresample = nresample)
    )
    
    return(structure(list(
      statistic  = c("Cochran's Q" = as.numeric(coin::statistic(res))),
      parameter  = c(df = NA_integer_),
      p.value    = as.numeric(coin::pvalue(res)),
      method     = paste0("Cochran's Q test (approximate, B = ", nresample, ")"),
      data.name  = DNAME
    ), class = "htest"))
  }
  
  ## -------------------------------------------------------------------
  ## Asymptotic
  ## -------------------------------------------------------------------
  k     <- nlevels(groups)
  y_mat <- matrix(unlist(split(y, blocks)), ncol = k, byrow = TRUE)
  y_mat <- y_mat[complete.cases(y_mat), , drop = FALSE]
  
  .cochranAsymptotic(y_mat, DNAME)
  
}



#' @rdname cochranQTest
#' @export
cochranQTest.formula <- function(formula,
                                 data,
                                 subset,
                                 na.action = na.pass,
                                 method    = c("asymptotic", "approximate"),
                                 nresample = 1e4,
                                 ...) {
  
  if (missing(formula) || length(formula) != 3L)
    stop("'formula' missing or incorrect")
  
  args <- list(
    formula   = formula,
    na.action = na.action,
    allowed   = "n.sample.dependent"
  )
  
  if (!missing(data))
    args$data <- data
  
  if (!missing(subset))
    args$subset <- substitute(subset)
  
  d <- do.call(bedrock::resolveFormula, args)
  
  res <- cochranQTest.default(
    y         = d$response,
    groups    = d$group,
    blocks    = d$block,
    method    = method,
    nresample = nresample,
    ...
  )
  
  res$data.name <- d$data.name
  res
}




# == internal helper functions ============================================


.cochranAsymptotic <- function(x, DNAME) {
  
  n <- nrow(x)
  k <- ncol(x)
  
  if (k < 2) stop("need at least 2 groups")
  if (n < 2) stop("need at least 2 blocks")
  
  Sj <- colSums(x)
  Ri <- rowSums(x)
  T  <- sum(x)
  
  denom <- k*T - sum(Ri^2)
  
  if (denom == 0) {
    # alle Blöcke haben identische Antworten -> Q = 0
    return(structure(list(
      statistic = c("Cochran's Q" = 0),
      parameter = c(df = k - 1L),
      p.value   = 1,
      method    = "Cochran's Q test (asymptotic)",
      data.name = DNAME
    ), class = "htest"))
  }
  
  
  Q <- (k - 1) * (k * sum(Sj^2) - T^2) / denom
  
  PVAL <- pchisq(Q, df = k - 1, lower.tail = FALSE)
  
  structure(list(
    statistic = c("Cochran's Q" = Q),
    parameter = c(df = k - 1),
    p.value = PVAL,
    method = "Cochran's Q test (asymptotic)",
    data.name = DNAME
  ), class = "htest")
}

