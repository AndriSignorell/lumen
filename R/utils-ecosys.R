

# internal utilities for the DescToolsX ecosystem
#    (to be used in every package, where needed)


#' Check and resolve verbose level
#'
#' @description
#' Resolves verbosity level using the following priority:
#' \itemize{
#'   \item function argument
#'   \item global option \code{DescTools.verbose}
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


#' ## ============================================================
#' ## Argument handling utilities (centralized, reusable)
#' ## ============================================================
#' 
#' #' @keywords internal
#' .extractBootArgs <- function(dots) {
#'   
#'   # NOTE: default type = "bca" here does not match what several CI
#'   # functions' own documentation states (e.g. meanCI/meanDiffCI/
#'   # quantileCI's docs say the default is "basic"). One of the two needs
#'   # to be reconciled -- either update this default to "basic", or update
#'   # those functions' @param \dots docs to say "bca".
#'   bedrock::extractArgs(
#'     dots,
#'     defaults = list(
#'       type     = "bca",
#'       R        = 999,
#'       parallel = "no",
#'       ncpus    = getOption("boot.ncpus", 1L)
#'     ),
#'     validate = function(x) {
#'       
#'       if (!is.numeric(x$R) || length(x$R) != 1 || x$R <= 0)
#'         stop("R must be a positive integer")
#'       
#'       if (!is.numeric(x$ncpus) || x$ncpus < 1)
#'         stop("ncpus must be >= 1")
#'       
#'       if (!x$parallel %in% c("no", "multicore", "snow"))
#'         stop("parallel must be 'no', 'multicore', or 'snow'")
#'     }
#'   )
#' }
