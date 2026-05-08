
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
#' @param x A contingency table (matrix or array of counts).
#' @param MARGIN An integer indicating the margin over which to compute
#'   the scores. Defaults to \code{1} (rows). Use \code{2} for columns.
#' @param method A character string specifying the scoring method.
#'   One of:
#'   \itemize{
#'     \item \code{"table"}: Uses numeric dimnames if available, otherwise
#'       assigns sequential integers.
#'     \item \code{"ranks"}: Mid-ranks based on cumulative frequencies.
#'     \item \code{"ridit"}: Ridit scores (ranks divided by total count).
#'     \item \code{"modridit"}: Modified ridit scores (ranks divided by
#'       total count + 1).
#'   }
#'
#' @details
#' For \code{method = "table"}, numeric dimension names are used as scores
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
#' \url{https://stat.ethz.ch/pipermail/r-help/2005-July/076371.html}
#'
#' @seealso \code{\link{cochranArmitageTest}}, \code{\link{cor}}
#' 


#' @family table.utils
#' @concept table-manipulation
#' @concept nonparametric
#' @concept descriptive-statistics
#'
#'
#' @export
scores <- function(x, MARGIN=1, 
                   method=c("table", "ranks", "ridit", "modridit")) { 
  
  # used by cochranArmitageTest, pearsonCor, spearmanCor
  
  
  # original by Eric Lecoutre
  # https://stat.ethz.ch/pipermail/r-help/2005-July/076371.html
  
  if (method == "table"){
    
    if (is.null(dimnames(x)) || 
        any(is.na(suppressWarnings(as.numeric(dimnames(x)[[MARGIN]]))))) {
      res <- 1:dim(x)[MARGIN]
    } else {
      res <- (as.numeric(dimnames(x)[[MARGIN]]))
    }
    
  } else	{
    ### method is a rank one
    Ndim <- dim(x)[MARGIN]
    OTHERMARGIN <- 3 - MARGIN
    
    ranks <- c(0, (cumsum(apply(x, MARGIN, sum))))[1:Ndim] + 
      (apply(x, MARGIN, sum)+1) /2 
    
    if (method == "ranks") res <- ranks
    if (method == "ridit") res <- ranks/(sum(x))
    if (method == "modridit") res <- ranks/(sum(x)+1)
  }
  
  return(res)
  
}

