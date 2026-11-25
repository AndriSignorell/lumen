
#' Compute Scores for Ordinal Contingency Tables
#'
#' Computes score values for the levels of a contingency table margin.
#' These scores are used in several statistical procedures such as the
#' Cochran-Armitage test and correlation measures for ordinal data.
#'
#' The function supports different scoring methods, including simple
#' table-based scores, ranks, and ridit-type transformations.
#'
#' @param x a contingency table (matrix or array of counts).
#' @param margin an integer indicating the margin over which to compute
#'   the scores. Defaults to `1` (rows). Use `2` for columns.
#' @param method a character string specifying the scoring method.
#'   One of:
#'   \itemize{
#'     \item `"table"`: uses numeric dimnames if available, otherwise
#'       assigns sequential integers.
#'     \item `"ranks"`: mid-ranks based on cumulative frequencies.
#'     \item `"ridit"`: ridit scores (ranks divided by total count).
#'     \item `"mod-ridit"`: modified ridit scores (ranks divided by
#'       total count + 1).
#'   }
#'
#' @details
#' For `method = "table"`, numeric dimension names are used as scores
#' if available. Otherwise, consecutive integers starting from 1 are assigned.
#'
#' For rank-based methods, scores are computed as midpoints of cumulative
#' frequencies along the selected margin.
#'
#' Ridit and modified ridit scores are normalized versions of these ranks.
#'
#' @return A numeric vector of scores corresponding to the levels of the
#'   selected margin.
#'
#' @references
#' Lecoutre, E. (2005). R-help mailing list discussion.
#' <https://stat.ethz.ch/pipermail/r-help/2005-July/076371.html>
#'
#' @seealso [cochranArmitageTest()], [cor()]
#'
#' @family scores
#' @concept transformation
#' @concept ordinal
#'
#' @examples
#' # dose-response table: numeric dimnames are used as scores
#' tab <- matrix(c(10, 9, 10, 7,
#'                  0, 1,  0, 3), byrow = TRUE, nrow = 2,
#'               dimnames = list(resp = 0:1, dose = c(0, 1, 2, 5)))
#' scores(tab, margin = 2)                    # 0 1 2 5
#'
#' # rank-based scores use the column totals (10 each)
#' scores(tab, margin = 2, method = "ranks")  # 5.5 15.5 25.5 35.5
#' scores(tab, margin = 2, method = "ridit")
#'
#' # non-numeric dimnames fall back to 1, 2, ...
#' dimnames(tab)$dose <- c("none", "low", "mid", "high")
#' scores(tab, margin = 2)
#'
#' @export
scores <- function(x, margin = 1,
                   method = c("table", "ranks", "ridit", "mod-ridit")) {

  # used by cochranArmitageTest, pearsonCor, spearmanCor
  # original by Eric Lecoutre
  # https://stat.ethz.ch/pipermail/r-help/2005-July/076371.html

  method <- match.arg(method)

  if (method == "table") {
    lev <- dimnames(x)[[margin]]
    num <- suppressWarnings(as.numeric(lev))
    if (is.null(lev) || anyNA(num))
      return(seq_len(dim(x)[margin]))
    return(num)
  }

  # rank-based methods: mid-ranks of the cumulative margin totals
  tot   <- apply(x, margin, sum)
  ranks <- cumsum(tot) - tot + (tot + 1) / 2

  switch(method,
         "ranks"     = ranks,
         "ridit"     = ranks / sum(x),
         "mod-ridit" = ranks / (sum(x) + 1))
}
