

## == internal helper functions ==========================================
## Note: take care for file name starting aaa to ensure first load!



utils::globalVariables(c(
  "zTest", "runsTest", "varTest", 
  "pageTest", "jonckheereTerpstraTest", "cochranQTest", "siegelTukeyTest"
))




#' @keywords internal
.replace_text_calls <- function(expr, old = "oldFun", new = "newFun") {
  
  if (is.call(expr)) {
    
    # case 1: simple symbol, e.g. t.test(...)
    if (is.symbol(expr[[1]]) && as.character(expr[[1]]) == old) {
      expr[[1]] <- as.name(new)
    }
    
    # case 2: namespace qualified, e.g. stats::t.test(...)
    if (is.call(expr[[1]]) &&
        identical(expr[[1]][[1]], as.name("::")) &&
        as.character(expr[[1]][[3]]) == old) {
      expr[[1]][[3]] <- as.name(new)
    }
    
    # case 3: do.call("t.test", ...)
    if (identical(expr[[1]], as.name("do.call")) &&
        is.character(expr[[2]]) && expr[[2]] == old) {
      expr[[2]] <- new
    }
    
    # recursively over all arguments (not expr[[1]], this has been handled before)
    expr[-1] <- lapply(as.list(expr[-1]), .replace_text_calls, old = old, new = new)
  }
  
  expr
}




#' @keywords internal
.asBinary <- function(x, ref=NULL, warn = TRUE) {
  
  # remove names (safer downstream)
  x <- unname(x)
  
  # -------------------------
  # logical
  # -------------------------
  if (is.logical(x)) {
    return(as.numeric(x))
  }
  
  # -------------------------
  # numeric
  # -------------------------
  if (is.numeric(x)) {
    if (!all(x %in% c(0,1, NA)))
      stop("'x' must be binary (0/1)")
    return(x)
  }
  
  # -------------------------
  # factor
  # -------------------------
  if (is.factor(x)) {
    
    lev <- levels(x)
    
    if (length(lev) != 2)
      stop("'x' must have exactly 2 levels")
    
    # explicit reference level
    if (!is.null(ref)) {
      
      if (!ref %in% lev)
        stop("'ref' must be one of the factor levels")
      
      return(as.numeric(x == ref))
    }
    
    # default behavior
    if (warn)
      warning(gettextf("coercing factor to binary (0/1): using '%s' as '1'", lev[2]))
    
    return(as.numeric(x == lev[2]))
  }
  
  
  # -------------------------
  # character
  # -------------------------
  if (is.character(x)) {
    u <- unique(x)
    u <- u[!is.na(u)]
    
    if (length(u) != 2)
      stop("'x' must have exactly 2 unique values")
    
    # explicit reference level
    if (!is.null(ref)) {
      
      if (!ref %in% u)
        stop("'ref' must be one of the unique values")
      
      return(as.numeric(x == ref))
    }
    
    # default behavior
    if (warn)
      warning(gettextf("coercing factor to binary (0/1): using '%s' as '1'", u[2]))
    
    return(as.numeric(x == u[2]))
    
  }
  
  # -------------------------
  # fallback
  # -------------------------
  stop("unsupported type for 'x'")
}

