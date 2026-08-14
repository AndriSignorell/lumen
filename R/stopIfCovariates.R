

#' Stop on covariates in the model
#'
#' The pairwise procedures take their group means from model.tables(), which
#' projects sequentially and therefore reports means unadjusted for any
#' covariate; they would also depend on the order of the terms in the formula.
#' The standard errors would additionally miss the term for the difference of
#' the covariate means.
#'
#' @keywords internal
#' @noRd
.stopIfCovariates <- function(x) {
  
  tt <- stats::terms(x)
  cls <- attr(tt, "dataClasses")
  
  if (is.null(cls)) {
    return(invisible(NULL))
  }
  
  resp <- attr(tt, "response")
  
  if (resp > 0L) {
    cls <- cls[-resp]
  }
  
  # nmatrix.k covers poly() and matrix covariates
  isCov <- grepl("^(numeric|nmatrix)", cls)
  
  if (any(isCov)) {
    stop(gettextf(
      "the model contains the covariate(s) %s, the group means would not be adjusted for them; use emmeans::emmeans() for post hoc comparisons after ancova",
      paste(sQuote(names(cls)[isCov]), collapse = ", ")
    ), call. = FALSE)
  }
  
  invisible(NULL)
}

