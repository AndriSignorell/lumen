
# == internal helper functions ============================================
#
# Shared by every test that accepts an external vector alongside a model
# frame (orderBy, blocks, weights). They belong to lumen, not to bedrock:
# resolveFormula() hands out 'rows', these two turn it into something a
# test can use.


#' Positions of the retained rows in the original data
#'
#' Counterpart of the `rows` component of [resolveFormula()], for objects
#' that carry a model frame of their own instead of one built here, such as
#' a fitted `"lm"`. The row names are the only representation surviving both
#' `subset` and `na.action`; they are matched against `data` when it is
#' available and only read as row numbers as a last resort, so that
#' non-default row names give `NULL`, never a wrong alignment.
#'
#' @param rn character, the row names of the model frame or design matrix.
#' @param data the data the model was fitted from, if available.
#'
#' @return an integer vector of positions, or `NULL` if they cannot be
#'   determined.
#'
#' @noRd
.rowsFromNames <- function(rn, data = list()) {

  if (is.null(rn))
    return(NULL)

  idx <- if (!is.null(rownames(data)))
    match(rn, rownames(data))
  else
    suppressWarnings(as.integer(rn))

  if (anyNA(idx)) NULL else idx
}


#' Permutation index for an external ordering variable
#'
#' Resolves `orderBy` against the data, aligns it with the model frame and
#' returns the permutation. The caller applies it to whatever it holds, so
#' the helper stays independent of the number and shape of the objects to be
#' reordered:
#'
#' \preformatted{
#' ord <- .orderIndex(orderBy, nrow(X), data = data, rows = r$rows)
#' if (!is.null(ord)) {
#'   X <- X[ord, , drop = FALSE]
#'   y <- y[ord]
#' }
#' }
#'
#' An ordering variable is written against the original data, while the
#' model frame has been filtered by `subset` and by `na.action`. Whenever
#' the lengths disagree, `rows` is used to reduce the ordering variable to
#' exactly the retained observations; without usable `rows` this is an error
#' rather than a silently wrong ordering.
#'
#' Missing values in the ordering variable are left to [order()] and end up
#' last.
#'
#' @param orderBy a vector, a data frame, a one-sided formula, or `NULL`.
#'   Several columns or terms are used as successive ordering keys. Each key
#'   must hold one value per observation.
#' @param n integer, the number of rows to be ordered.
#' @param data an optional data frame the formula is evaluated in.
#' @param rows integer, the positions of the retained rows in the original
#'   data, as returned by [bedrock::resolveFormula()] or
#'   `NULL`.
#'
#' @return an integer permutation of length `n`, or `NULL` if `orderBy` is
#'   `NULL`.
#'
#' @noRd
.orderIndex <- function(orderBy, n, data = list(), rows = NULL) {

  if (is.null(orderBy))
    return(NULL)

  # a list of ordering keys, whatever the input was
  keys <- if (inherits(orderBy, "formula"))
    unname(as.list(model.frame(orderBy, data = data, na.action = na.pass)))
  else if (is.data.frame(orderBy))
    unname(as.list(orderBy))
  else
    list(orderBy)

  if (!length(keys) || !length(keys[[1L]]))
    stop("'orderBy' is empty", call. = FALSE)

  # A key must hold one value per observation. Both sources above guarantee
  # equal key lengths, so checking the first one is enough - but not that a
  # key is a plain vector: a term such as ~ poly(t, 2) yields a single
  # matrix column, which would be flattened into 2n values and produce an
  # ordering longer than the model frame.
  if (any(vapply(keys, function(z) !is.null(dim(z)), NA)))
    stop("'orderBy' must have one value per observation", call. = FALSE)

  len <- length(keys[[1L]])

  if (len != n) {

    if (is.null(rows) || length(rows) != n || max(rows) > len)
      stop(gettextf(
        "'orderBy' cannot be aligned with the %d rows of the model frame", n),
        call. = FALSE)

    keys <- lapply(keys, `[`, rows)
  }

  # do.call() is safe here: keys holds plain atomic vectors, never an
  # unevaluated expression
  do.call(order, keys)
}
