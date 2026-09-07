#' Print Method for rankTest Objects
#'
#' Prints pairwise comparison results produced by [dunnTest()],
#' [conoverTest()], or [nemenyiTest()].
#'
#' @param x an object of class `"rankTest"`.
#' @param digits number of significant digits used for printing numeric values.
#'   Passed to [print.data.frame()]. Defaults to
#'   `getOption("digits", 3)`.
#' @param \dots further arguments passed to [print.data.frame()] or
#'   [print.default()].
#'
#' @return `x`, invisibly.
#'
#' @seealso [dunnTest()], [conoverTest()],
#'   [nemenyiTest()]
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
