#' Print Method for rankTest Objects
#'
#' Prints pairwise comparison results produced by \code{\link{dunnTest}},
#' \code{\link{conoverTest}}, or \code{\link{nemenyiTest}}.
#'
#' @param x an object of class \code{"rankTest"}.
#' @param digits number of significant digits used for printing numeric values.
#'   Passed to \code{\link{print.data.frame}}. Defaults to
#'   \code{getOption("digits", 3)}.
#' @param \dots further arguments passed to \code{\link{print.data.frame}} or
#'   \code{\link{print.default}}.
#'
#' @return \code{x}, invisibly.
#'
#' @seealso \code{\link{dunnTest}}, \code{\link{conoverTest}},
#'   \code{\link{nemenyiTest}}
#'

#' @export
print.rankTest <- function(
    x,
    digits = getOption("digits", 3),
    ...
) {

  cat("\n", attr(x, "main"), "\n\n")

  if (attr(x, "output") == "list") {

    xx <- data.frame(x$res)

    xx$" " <- fm(xx$pval, fmt = "*")

    xx$pval <- format.pval(
      xx$pval,
      digits = 2,
      nsmall = 4
    )

    print.data.frame(xx, digits = digits, ...)
    .printSignifCodes()

  } else {

    xx <- x$res

    xx[] <- format.pval(
      xx,
      digits = 2,
      na.form = "-"
    )

    print(xx, digits = digits, quote = FALSE, ...)
  }

  cat("\n")

  invisible(x)
}
