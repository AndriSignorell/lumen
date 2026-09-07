
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
