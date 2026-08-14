
#' Significantly Different Levels per Group
#'
#' Extracts, for every factor level, the set of levels it differs
#' significantly from, based on the pairwise comparisons of a post hoc test.
#'
#' The result answers the question a reader of a post hoc table usually has:
#' not which of the pairs are significant, but which levels a given level is
#' distinguishable from, and in which direction.
#'
#' The labels are numbers by default. Letters are available, but carry the
#' opposite meaning of a compact letter display, where levels sharing a letter
#' are those that could \emph{not} be distinguished; here a label marks a
#' significant difference.
#'
#' @param x an object of class \code{PostHocTest} as returned by
#'   \code{\link{postHocTest}}, or of class \code{pairwise.htest} as returned
#'   by \code{\link[stats]{pairwise.t.test}} and friends
#' @param alpha the significance level; defaults to \code{1 - conf.level} of
#'   the object, or to 0.05 where no confidence level is stored
#' @param direction logical; if \code{TRUE}, a sign is appended to every
#'   label, giving the direction of the row level's mean relative to the
#'   listed level. Requires the differences to be present in the object, which
#'   is not the case for the p-value branch and for \code{pairwise.htest}
#'   objects.
#' @param labels either a keyword defining how the levels are abbreviated, one
#'   of \code{"numbers"} (the default), \code{"letters"}, \code{"LETTERS"},
#'   \code{"abbreviate"} or \code{"names"}, or a character vector of labels
#'   with one entry per level
#' @param minlength the minimum length of the abbreviations, used for
#'   \code{labels = "abbreviate"} only
#' @param sep the separator between the labels
#' @param \dots further arguments, not used so far
#'
#' @return a list with one data frame per term, each holding the columns
#'   \code{label} (the label of the level itself) and \code{diff} (the labels
#'   of the levels it differs significantly from), with the levels as row
#'   names. The class is \code{"signifDiff"}, the significance level is kept
#'   in the attribute \code{alpha}.
#'
#' @examples
#' r.aov <- aov(breaks ~ tension, data = warpbreaks)
#' res <- postHocTest(r.aov, method = "hsd")
#'
#' signifDiff(res)
#'
#' # stricter level, without recomputing the test
#' signifDiff(res, alpha = 0.01)
#'
#' # readable labels instead of numbers
#' signifDiff(res, labels = "abbreviate", minlength = 4)
#'
#' # the p-value branch carries no differences, hence no signs
#' signifDiff(postHocTest(r.aov, method = "hsd", conf.level = NA))
#'
#' signifDiff(pairwise.t.test(warpbreaks$breaks, warpbreaks$tension))
#'
#' @seealso \code{\link{postHocTest}}, \code{\link[stats]{pairwise.t.test}}
#'
#' @family test.posthoc
#' @concept multiple-comparisons
#' @concept post-hoc
#'
#' @export
signifDiff <- function(x, ...) {
  UseMethod("signifDiff")
}




#' @rdname signifDiff
#' @export
signifDiff.PostHocTest <- function(x,
                                   alpha = NULL,
                                   direction = TRUE,
                                   labels = "numbers",
                                   minlength = 3L,
                                   sep = ", ",
                                   ...) {

  if (is.null(alpha)) {
    conf.level <- attr(x, "conf.level")

    alpha <- if (is.null(conf.level) || is.na(conf.level)) {
      0.05
    } else {
      1 - conf.level
    }
  }

  mats <- lapply(unclass(x), .pairMatrices)

  # only complain where the direction was asked for explicitly, the p-value
  # branch never holds differences
  if (!missing(direction) && isTRUE(direction) &&
      any(vapply(mats, function(m) is.null(m$d), NA))) {
    warning("the object holds no differences, the direction is dropped",
            call. = FALSE)
  }

  res <- lapply(mats, .diffTable, alpha = alpha, direction = direction,
                labels = labels, minlength = minlength, sep = sep)

  structure(
    res,
    alpha = alpha,
    signed = isTRUE(attr(res[[1L]], "signed")),
    method.str = attr(x, "method.str"),
    class = "signifDiff"
  )
}




#' @rdname signifDiff
#' @export
signifDiff.pairwise.htest <- function(x,
                                      alpha = 0.05,
                                      direction = FALSE,
                                      labels = "numbers",
                                      minlength = 3L,
                                      sep = ", ",
                                      ...) {

  if (!missing(direction) && isTRUE(direction)) {
    warning("a pairwise.htest holds no differences, the direction is dropped",
            call. = FALSE)
  }

  res <- list(
    .diffTable(.pairMatrices(x$p.value), alpha = alpha, direction = FALSE,
               labels = labels, minlength = minlength, sep = sep)
  )

  names(res) <- x$data.name

  structure(
    res,
    alpha = alpha,
    signed = FALSE,
    method.str = gettextf("\n  Pairwise comparisons : %s \n", x$method),
    class = "signifDiff"
  )
}




#' @rdname signifDiff
#' @param legend logical; if \code{TRUE}, the meaning of the signs is printed
#'   as a footer, separated by a rule
#' @export
print.signifDiff <- function(x, legend = TRUE, ...) {

  cat(attr(x, "method.str"))

  cat(gettextf("    levels a level differs from, at alpha = %s\n",
               format(attr(x, "alpha"))))

  for (nm in names(x)) {
    cat("\n$", nm, "\n", sep = "")
    print.data.frame(x[[nm]], right = FALSE, ...)
  }

  # constant over all terms, hence printed once, as the significance codes are
  if (legend && isTRUE(attr(x, "signed"))) {
    cat("\n---\n")
    cat("Sign codes:  '+' mean of the row level is higher than the listed one, '-' lower\n")
  }

  cat("\n")

  invisible(x)
}




# == internal helper functions ================================================


#' Symmetric matrices of the pairwise p-values and differences
#'
#' Copes with the two shapes a post hoc result comes in: the four column array
#' of the confidence interval branch, and the truncated triangular matrix of
#' the p-value branch, which pairwise.htest objects share. The differences are
#' available in the first case only, then \code{d[i, j]} is the mean of level
#' \code{i} minus the mean of level \code{j}.
#'
#' @keywords internal
#' @noRd
.pairMatrices <- function(z) {

  if (is.matrix(z) && "pval" %in% colnames(z)) {

    lvls <- attr(z, "levels")

    if (is.null(lvls)) {
      lvls <- .deduceLevels(rownames(z))
    }

    k <- length(lvls)

    if (k * (k - 1) / 2 != nrow(z)) {
      stop("the number of comparisons does not match the number of levels",
           call. = FALSE)
    }

    empty <- matrix(NA_real_, nrow = k, ncol = k, dimnames = list(lvls, lvls))

    # the pairs are stored in the column-wise order of lower.tri(), so the
    # index pairs can be reconstructed without parsing any row name
    idx <- which(lower.tri(empty), arr.ind = TRUE)
    rev <- idx[, c(2L, 1L), drop = FALSE]

    p <- empty
    p[idx] <- z[, "pval"]
    p[rev] <- z[, "pval"]

    d <- empty
    d[idx] <- z[, "diff"]
    d[rev] <- -z[, "diff"]

    return(list(p = p, d = d))
  }

  # truncated triangle: rows hold the levels 2..k, columns the levels 1..k-1
  lvls <- c(colnames(z), rownames(z)[nrow(z)])

  k <- length(lvls)

  p <- matrix(NA_real_, nrow = k, ncol = k, dimnames = list(lvls, lvls))

  p[2L:k, 1L:(k - 1L)] <- as.matrix(z)

  # the branch stores the lower triangle, mirror it
  p[upper.tri(p)] <- t(p)[upper.tri(p)]

  list(p = p, d = NULL)
}




#' Recover the level names from the pair labels
#'
#' Fallback for objects carrying no levels attribute. Exploits the order in
#' which the pairs were built: the second parts run level 1 (k-1 times),
#' level 2 (k-2 times) and so on, the last first part is the final level.
#' Fails on level names containing the separator, which is exactly why
#' postHocTest() should store the levels.
#'
#' @keywords internal
#' @noRd
.deduceLevels <- function(nms) {

  parts <- strsplit(nms, "-", fixed = TRUE)

  if (any(lengths(parts) != 2L)) {
    stop("the level names cannot be recovered from the pair labels, ",
         "recompute the test with the current version of postHocTest()",
         call. = FALSE)
  }

  first <- vapply(parts, `[`, "", 1L)
  second <- vapply(parts, `[`, "", 2L)

  c(unique(second), first[length(first)])
}




#' Turn the pairwise matrices into the level-wise listing
#'
#' @keywords internal
#' @noRd
.diffTable <- function(m, alpha, direction, labels, minlength, sep) {

  p <- m$p

  lvls <- rownames(p)

  lab <- .makeLabels(lvls, labels = labels, minlength = minlength)

  isSignif <- !is.na(p) & p < alpha

  diag(isSignif) <- FALSE

  signed <- isTRUE(direction) && !is.null(m$d)

  cells <- if (signed) {
    matrix(
      paste0(rep(lab, each = length(lab)), ifelse(m$d > 0, "+", "-")),
      nrow = length(lab)
    )
  } else {
    matrix(rep(lab, each = length(lab)), nrow = length(lab))
  }

  res <- data.frame(
    label = lab,
    diff = vapply(
      seq_along(lvls),
      function(i) paste(cells[i, isSignif[i, ]], collapse = sep),
      ""
    ),
    row.names = lvls,
    stringsAsFactors = FALSE
  )

  attr(res, "signed") <- signed

  res
}




#' @keywords internal
#' @noRd
.makeLabels <- function(lvls, labels, minlength = 3L) {

  k <- length(lvls)

  if (length(labels) > 1L) {

    if (length(labels) != k) {
      stop(gettextf("'labels' must have one entry per level (%d)", k),
           call. = FALSE)
    }

    return(as.character(labels))
  }

  labels <- match.arg(
    labels,
    c("numbers", "letters", "LETTERS", "abbreviate", "names")
  )

  switch(
    labels,
    numbers = as.character(seq_len(k)),
    letters = .letterSeq(k, letters),
    LETTERS = .letterSeq(k, LETTERS),
    abbreviate = unname(abbreviate(lvls, minlength = minlength)),
    names = lvls
  )
}




#' Label sequence carrying on with aa, ab, ... beyond the alphabet
#'
#' @keywords internal
#' @noRd
.letterSeq <- function(k, alphabet = letters) {

  n <- length(alphabet)

  if (k <= n) {
    return(alphabet[seq_len(k)])
  }

  res <- alphabet

  while (length(res) < k) {
    res <- c(res, paste0(rep(alphabet, each = n), alphabet))
  }

  res[seq_len(k)]
}
