
# == validation helpers for the distribution and moment functions ============
#
# Three helpers cover everything the d-p-q-r and m functions need to check.
# They are deliberately terse: the logical flags (log, lower.tail, log.p) are
# not validated, matching the base R distribution functions.


# A single value within an admissible range.
#
# The moment functions are scalar throughout (they branch on the parameters
# to decide whether a moment exists), and a few distribution functions
# require a scalar too. The numeric, finite and whole-number checks are
# delegated to isNumeric(), which also supplies the tolerance for the latter.
#
# @param x the value to check.
# @param lower,upper bounds of the admissible range.
# @param strictLower,strictUpper whether the respective bound is excluded.
# @param integerValued whether the value must be a whole number.
# @param name the argument name to report; taken from the call by default.
#
# @return `x`, invisibly. Called for the side effect of raising an error.
#
# @noRd
.assertScalar <- function(x, lower = -Inf, upper = Inf,
                          strictLower = FALSE, strictUpper = FALSE,
                          integerValued = FALSE,
                          name = deparse1(substitute(x))) {

  if (length(x) != 1L || !isNumeric(x))
    stop(gettextf("'%s' must be a single finite numeric value", name),
         call. = FALSE)

  if (integerValued && !isNumeric(x, isIntegerValued = TRUE))
    stop(gettextf("'%s' must be a whole number", name), call. = FALSE)

  if ((if (strictLower) x <= lower else x < lower) ||
      (if (strictUpper) x >= upper else x > upper))
    stop(gettextf("'%s' must lie in %s%s, %s%s", name,
                  if (strictLower) "(" else "[",
                  format(lower), format(upper),
                  if (strictUpper) ")" else "]"),
         call. = FALSE)

  invisible(x)
}


# A parameter that must be positive throughout, possibly a vector.
#
# Missing values are passed on rather than rejected, so that they propagate
# to the result as they do in the base R distribution functions.
#
# @noRd
.assertPositive <- function(x, name = deparse1(substitute(x))) {

  if (!is.numeric(x))
    stop(gettextf("'%s' must be numeric", name), call. = FALSE)

  if (any(x <= 0, na.rm = TRUE))
    stop(gettextf("'%s' must be positive", name), call. = FALSE)

  invisible(x)
}


# Bring the 'p' of a quantile function onto the lower-tail probability scale.
#
# The boundary values 0 and 1 are admissible and yield the end points of the
# support; values outside [0, 1] give NaN with a warning, as in base R.
#
# @noRd
.qProb <- function(p, lower.tail = TRUE, log.p = FALSE) {

  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p

  if (any(bad <- !is.na(p) & (p < 0 | p > 1))) {
    warning("NaNs produced")
    p[bad] <- NaN
  }

  p
}
