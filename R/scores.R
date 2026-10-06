
#' Compute Scores for Ordinal Contingency Tables
#'
#' A utility function computing score transformations of raw data, 
#' including normal scores, exponential scores, and Savage scores, 
#' typically used as a preprocessing step for nonparametric tests.
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
#' @export
scores <- function(x, margin=1, 
                   method=c("table", "ranks", "ridit", "mod-ridit")) { 
  
  # used by cochranArmitageTest, pearsonCor, spearmanCor
  
  method <- match.arg(method)
  
  # original by Eric Lecoutre
  # https://stat.ethz.ch/pipermail/r-help/2005-July/076371.html
  
  if (method == "table"){
    
    if (is.null(dimnames(x)) || 
        any(is.na(suppressWarnings(as.numeric(dimnames(x)[[margin]]))))) {
      res <- 1:dim(x)[margin]
    } else {
      res <- (as.numeric(dimnames(x)[[margin]]))
    }
    
  } else	{
    ### method is a rank one
    Ndim <- dim(x)[margin]
    OTHERMARGIN <- 3 - margin
    
    ranks <- c(0, (cumsum(apply(x, margin, sum))))[1:Ndim] + 
      (apply(x, margin, sum)+1) /2 
    
    if (method == "ranks") res <- ranks
    if (method == "ridit") res <- ranks/(sum(x))
    if (method == "mod-ridit") res <- ranks/(sum(x)+1)
  }
  
  return(res)
  
}

