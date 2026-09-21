
#' Two-Sided Level Behind a One-Sided Bound
#'
#' Internal helper shared by the `binom*CI()` family. A one-sided bound at
#' level \eqn{\gamma} is the corresponding end of the two-sided interval at
#' level \eqn{2\gamma - 1} (design_rules 4.1), so the tail probability the
#' methods have to work with is \eqn{2(1 - \gamma)} rather than
#' \eqn{1 - \gamma}.
#'
#' @param conf.level the confidence level as requested by the user.
#' @param sides one of `"two.sided"`, `"left"` or `"right"`, already matched.
#'
#' @return a single numeric giving the alpha the interval methods use.
#'
#' @keywords internal
.sidesAlpha <- function(conf.level, sides)
  1 - if (sides == "two.sided") conf.level else 2 * conf.level - 1



#' Validate the level of a possibly one-sided interval
#'
#' Internal helper for the CI functions. A one-sided bound at level
#' \eqn{\gamma} is computed as the end of the two-sided interval at level
#' \eqn{2\gamma - 1} (design_rules 4.1), which needs \eqn{\gamma > 0.5};
#' below that the transformed level is zero or negative and the quantile
#' functions returned NaN bounds without complaint.
#'
#' @param conf.level the confidence level as requested by the user.
#' @param sides one of `"two.sided"`, `"left"` or `"right"`, already matched.
#'
#' @return `conf.level`, invisibly.
#'
#' @keywords internal
#' @noRd
.checkSidedLevel <- function(conf.level, sides) {

  if (!is.numeric(conf.level) || length(conf.level) != 1L ||
      !is.finite(conf.level) || conf.level <= 0 || conf.level >= 1)
    stop("'conf.level' must be a single number between 0 and 1",
         call. = FALSE)

  if (sides != "two.sided" && conf.level <= 0.5)
    stop(gettextf("a one-sided interval needs 'conf.level' above 0.5, not %g",
                  conf.level), call. = FALSE, domain = NA)

  invisible(conf.level)
}
