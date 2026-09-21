

# internal utilities for the DescToolsX ecosystem
#    (to be used in every package, where needed)


#' Check and resolve verbose level
#'
#' @description
#' Resolves verbosity level using the following priority:
#' \itemize{
#'   \item function argument
#'   \item global option `DescTools.verbose`
#'   \item default (2)
#' }
#'
#' @param verbose Optional integer (1-3).
#'
#' @return Integer in \{1, 2, 3\}.
#' @keywords internal
.checkVerbose <- function(verbose = NULL){
  
  # resolve: arg > option > default
  verbose <- if(!is.null(verbose)) {
    verbose
  } else {
    getOption("DescTools.verbose", 2L)
  }
  
  # validation
  if(length(verbose) != 1 || !is.numeric(verbose) || !(verbose %in% 1:3)){
    stop("verbose must be a single integer: 1 (minimal), 2 (standard), or 3 (detailed).")
  }
  
  as.integer(verbose)
}



#' @keywords internal
.printSignifCodes <- function() {
  cat(
    "---\n",
    "Signif. codes:  ",
    "0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n",
    sep = ""
  )
}


