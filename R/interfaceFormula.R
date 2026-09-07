
#' Formula Interface - Common Arguments
#'
#' Common formula-based interface shared by multiple functions.
#'
#' @name Formulas
#' 
#' @param formula a formula of the form `lhs ~ rhs`, where `lhs`
#'   gives the response values and `rhs` the corresponding groups
#'   or explanatory variables.
#'
#' @param data an optional matrix or data frame (or similar; see
#'   [model.frame()]) containing the variables in the
#'   formula. By default the variables are taken from
#'   `environment(formula)`.
#'
#' @param subset an optional vector specifying a subset of observations
#'   to be used in the analysis.
#'
#' @param na.action a function which indicates what should happen when
#'   the data contain `NA`s. Defaults to
#'   `getOption("na.action")`.
#'
#' @details
#' Formula interfaces are evaluated using [model.frame()],
#' following standard R conventions.
#' The left-hand side of the formula must contain the response variable.
#' The right-hand side typically specifies a grouping or explanatory variable.
#' Only formulas with a single response and at least one explanatory
#' variable are supported.
#'
#' See also:
#' \itemize{
#'   \item [formula()]
#'   \item [model.frame()]
#'   \item [terms()]
#' }
#'
#' @keywords internal
#' 
#' @family topic.systemTools
#' @concept formula interface
#' @concept model specification
#' @concept API
NULL
