
#' Post Hoc Tests for ANOVA
#'
#' Provides a unified interface for several parametric post hoc tests
#' following a significant ANOVA. The function computes pairwise
#' comparisons of group means with different methods for controlling
#' the family-wise error rate.
#'
#' \bold{Overview.}
#' Post hoc tests differ in how strongly they control type I error.
#' Conservative methods reduce false positives but have lower power,
#' while more liberal methods increase power at the cost of a higher
#' false positive rate.
#'
#' \bold{Available methods.}
#' \itemize{
#'   \item \strong{LSD (Fisher)}: No adjustment for multiple testing;
#'   highest power but inflated type I error. Mainly suitable for a small
#'   number of groups.
#'
#'   \item \strong{Bonferroni} (Dunn's (Bonferroni) t-test): Adjusts p-values 
#'   by the number of comparisons;
#'   simple and robust, but often overly conservative.
#'
#'   \item \strong{Tukey HSD}: Controls the family-wise error rate for all
#'   pairwise comparisons; widely used and generally recommended for balanced
#'   designs.
#'
#'   \item \strong{Newman-Keuls}: Stepwise procedure with more power than Tukey,
#'   but weaker error control; may inflate type I error.
#'
#'   \item \strong{Duncan}: Similar to Newman-Keuls but more liberal; provides
#'   higher power at the cost of increased false positives.
#'
#'   \item \strong{Scheffé}: Very conservative; suitable for both pairwise and
#'   complex (contrast-based) comparisons.
#' }
#'
#' \bold{Guidance.}
#' Tukey HSD is typically a good default for pairwise comparisons.
#' Bonferroni is useful when strict error control is required.
#' Scheffé is appropriate for more complex contrasts.
#'
#'
#' @return An object of class \code{PostHocTest}.
#'  
#' \bold{Tables} \verb{ } For tables pairwise chi-square test can be performed,
#' either without correction or with correction for multiple testing following
#' the logic in \code{\link{p.adjust}}. 
#' 
#' @name postHoc
#' @aliases postHocTest postHocTest.aov postHocTest.table postHocTest.matrix
#' print.postHocTest plot.postHocTest
#' @param x An object of class \code{aov}.
#' @param method one of \code{"hsd"}, \code{"bonf"}, \code{"lsd"},
#' \code{"scheffe"}, \code{"newmankeuls"}, defining the method for the pairwise
#' comparisons.\cr For the post hoc test of tables the methods of
#' \code{\link{p.adjust}} can be supplied. See the detail there. 
#' @param which a character vector listing terms in the fitted model for which
#' the intervals should be calculated. Defaults to all the terms. 
#' @param conf.level a numeric value between zero and one giving the
#' family-wise confidence level to use.  If this is set to NA, just a matrix
#' with the p-values will be returned. 
#' @param ordered a logical value indicating if the levels of the factor should
#' be ordered according to increasing average in the sample before taking
#' differences. If ordered is \code{TRUE} then the calculated differences in
#' the means will all be positive. The significant differences will be those
#' for which the lower end point is positive. \cr This argument will be ignored
#' if method is not either \code{hsd} or \code{newmankeuls}.
#' @param digits controls the number of fixed digits to print.
#' @param \dots further arguments, not used so far.
#'  
#' @return an object of type "postHocTest", which will either be \cr A) a list
#' of data.frames containing the mean difference, lower ci, upper ci and the
#' p-value, if a conf.level was defined (something else than NA) or \cr B) a
#' list of matrices with the p-values, if conf.level has been set to NA. 
#' 
#' @seealso \code{\link{TukeyHSD}}, \code{\link{aov}},
#' \code{\link{pairwise.t.test}}, \code{\link{scheffeTest}} 
#' 
#' @family topic.postHocTests
#' @concept multiple comparisons
#' @concept wrapper
#' @concept convenience function
#' 
#' @examples
#' 
#' postHocTest(aov(breaks ~ tension, data = warpbreaks), method = "lsd")
#' postHocTest(aov(breaks ~ tension, data = warpbreaks), method = "hsd")
#' postHocTest(aov(breaks ~ tension, data = warpbreaks), method = "scheffe")
#' 
#' r.aov <- aov(breaks ~ tension, data = warpbreaks)
#' 
#' # compare p-values:
#' round(cbind(
#'     lsd= postHocTest(r.aov, method="lsd")$tension[,"pval"]
#'   , bonf=postHocTest(r.aov, method="bonf")$tension[,"pval"]
#' ), 4)
#' 
#' # only p-values by setting conf.level to NA
#' postHocTest(aov(breaks ~ tension, data = warpbreaks), method = "hsd",
#'             conf.level=NA)
#' 


#' @rdname postHoc
#' @export
postHocTest <- function (x, ...)
  UseMethod("postHocTest")



#' @rdname postHoc
#' @export
postHocTest.aov <- function (x, which = NULL,
                             method=c("hsd","bonferroni","lsd","scheffe",
                                      "newmankeuls","duncan"),
                             conf.level = 0.95, ordered = FALSE, ...) {
  
  method <- match.arg(method)
  
  FUN_MAP <- list(
    bonferroni = .bonferroni,
    lsd = .lsd,
    hsd = .hsd,
    newmankeuls = .newmankeuls,
    duncan = .duncan,
    scheffe = .scheffe
  )
  
  mm <- model.tables(x, "means")
  if (is.null(mm$n))
    stop("no factors in the fitted model")
  
  tabs <- mm$tables[-1L]
  
  if(is.null(which)) which <- seq_along(tabs)
  tabs <- tabs[which]
  
  nn <- mm$n[names(tabs)]
  nn_na <- is.na(nn)
  
  if (all(nn_na))
    stop("'which' specified no factors")
  
  if (any(nn_na)) {
    warning("'which' specified some non-factors which will be dropped")
    tabs <- tabs[!nn_na]
    nn <- nn[!nn_na]
  }
  
  out <- setNames(vector("list", length(tabs)), names(tabs))
  MSE <- sum(x$residuals^2)/x$df.residual
  
  for (nm in names(tabs)) {
    
    tab <- tabs[[nm]]
    means <- as.vector(tab)
    
    nms <- if (length(dim(tab)) > 1L) {
      dn <- dimnames(tab)
      apply(do.call("expand.grid", dn), 1L, paste, collapse = ":")
    } else names(tab)
    
    n <- nn[[nm]]
    if (length(n) < length(means))
      n <- rep.int(n, length(means))
    
    if (method %in% c("hsd", "newmankeuls", "duncan") && isTRUE(ordered)) {
      ord <- order(means)
      means <- means[ord]
      n <- n[ord]
      if (!is.null(nms)) nms <- nms[ord]
    }
    
    center_mat <- outer(means, means, "-")
    keep <- lower.tri(center_mat)
    center <- center_mat[keep]
    
    k <- length(means)
    fun <- FUN_MAP[[method]]
    
    res <- if (method %in% c("newmankeuls", "duncan")) {
      fun(center = center, means = means, n = n,
          MSE = MSE, df = x$df.residual, conf.level = conf.level)
    } else {
      fun(center = center, n = n,
          MSE = MSE, df = x$df.residual,
          k = k, conf.level = conf.level)
    }
    
    width <- res$width
    pvals <- res$pvals
    method.str <- res$method.str
    
    if (!is.null(conf.level) && !is.na(conf.level)) {
      
      dnames <- list(NULL, c("diff", "lwr.ci", "upr.ci", "pval"))
      if (!is.null(nms))
        dnames[[1L]] <- outer(nms, nms, paste, sep = "-")[keep]
      
      out[[nm]] <- array(
        c(center, center - width, center + width, pvals),
        c(length(width), 4L),
        dnames
      )
      
    } else {
      
      mat <- matrix(NA, nrow = length(means), ncol = length(means))
      mat[lower.tri(mat)] <- pvals
      dimnames(mat) <- list(nms, nms)
      out[[nm]] <- mat[-1, -ncol(mat)]
    }
  }
  
  class(out) <- "PostHocTest"
  attr(out, "orig.call") <- x$call
  attr(out, "conf.level") <- conf.level
  attr(out, "ordered") <- ordered
  attr(out, "method") <- method.str
  attr(out, "method.str") <- gettextf(
    "\n  Posthoc multiple comparisons of means : %s \n",
    method.str
  )
  
  return(out)
}




#' @rdname postHoc
#' @export
postHocTest.matrix <- function(x, method = c("none","fdr","BH","BY","bonferroni","holm","hochberg","hommel"),
                               conf.level = 0.95, ...) {
  
  # http://support.sas.com/resources/papers/proceedings14/1544-2014.pdf
  
  # no conf.level supported so far
  conf.level  <- NA
  
  method <- match.arg(method)
  
  pvals <- pairApply(t(as.matrix(x)), 
                     FUN = function(y1, y2) 
                       chisq.test(cbind(y1,y2))$p.value, symmetric=TRUE)
  
  pvals[upper.tri(pvals, diag=TRUE)] <- NA
  
  if(method != "none")
    pvals[] <- p.adjust(pvals, method=method)
  
  #  pvals[] <- format.pval(pvals, digits = 2, na.form = "-")
  pvals <- pvals[-1, -ncol(pvals)]
  out <- list()
  out[[deparse(substitute(x))]] <- pvals
  
  class(out) <- c("PostHocTest")
  attr(out, "orig.call") <- "table"
  attr(out, "conf.level") <- conf.level
  attr(out, "ordered") <- FALSE
  attr(out, "method") <- method
  attr(out, "method.str") <- gettextf("\n  Posthoc multiple comparisons on chi-square test : %s \n", attr(out, "method"))
  
  return(out)
  
}


#' @rdname postHoc
#' @export
postHocTest.table <- function(x, method = c("none","fdr","BH","BY","bonferroni","holm","hochberg","hommel"),
                              conf.level = 0.95, ...) {
  class(x) <- "matrix"
  postHocTest(x, method=method, conf.level=conf.level, ...)
}




#' @rdname postHoc
#' @export
print.PostHocTest <- function(x, digits = getOption("digits", 3), ...) {
  
  method_str <- attr(x, "method.str")
  conf_level <- attr(x, "conf.level")
  ordered <- attr(x, "ordered")
  orig_call <- attr(x, "orig.call")
  
  cat(method_str)
  
  if (!is.null(conf_level) && !is.na(conf_level)) {
    cat("    ", format(100 * conf_level, digits = 2),
        "% family-wise confidence level\n", sep = "")
  }
  
  if (isTRUE(ordered)) {
    cat("    factor levels have been ordered\n")
  }
  
  if (!is.null(orig_call)) {
    cat("\nFit: ", deparse(orig_call, width.cutoff = 500L), "\n\n", sep = "")
  } else {
    cat("\n")
  }
  
  xx <- unclass(x)
  
  attributes(xx)[c("orig.call", "conf.level", "ordered", "method.str", "method")] <- NULL
  xx[["data.name"]] <- NULL
  
  if (!is.null(conf_level) && !is.na(conf_level)) {
    
    xx <- lapply(xx, function(xi) as.data.frame(xi))
    
    for (nm in names(xx)) {
      if ("pval" %in% names(xx[[nm]])) {
        xx[[nm]]$signif <- aurora::fm(xx[[nm]]$pval, fmt = "*")
        xx[[nm]]$pval <- aurora::fm(xx[[nm]]$pval, fmt = "p")
      }
    }
    
    print.default(xx, digits = digits, ...)
    
    cat("---\nSignif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
    
  } else {
    
    for (nm in names(xx)) {
      xx[[nm]][] <- aurora::fm(xx[[nm]], fmt = "p", na.form = "-")
    }
    
    print(xx, digits = digits, quote = FALSE, ...)
  }
  
  cat("\n")
  
  invisible(x)
}




# == internal helper functions ================================================


.bonferroni <- function(center, n, MSE, df, k, conf.level) {
  se <- sqrt(MSE * outer(1/n, 1/n, "+"))
  keep <- lower.tri(se)
  
  width <- qt(1 - (1 - conf.level)/(k * (k - 1)), df) * se[keep]
  est <- center / se[keep]
  pvals <- pmin(2 * pt(abs(est), df = df, lower.tail = FALSE) * (k*(k-1)/2), 1)
  
  list(width=width, est=est, pvals=pvals, method.str="Bonferroni")
}


.lsd <- function(center, n, MSE, df, k, conf.level) {
  se <- sqrt(MSE * outer(1/n, 1/n, "+"))
  keep <- lower.tri(se)
  
  width <- qt(1 - (1 - conf.level)/2, df) * se[keep]
  est <- center / se[keep]
  pvals <- 2 * pt(abs(est), df = df, lower.tail = FALSE)
  
  list(width=width, est=est, pvals=pvals, method.str="Fisher LSD")
}


.hsd <- function(center, n, MSE, df, k, conf.level) {
  se <- sqrt((MSE/2) * outer(1/n, 1/n, "+"))
  keep <- lower.tri(se)
  
  width <- qtukey(conf.level, k, df) * se[keep]
  est <- center / se[keep]
  pvals <- ptukey(abs(est), k, df, lower.tail = FALSE)
  
  list(width=width, est=est, pvals=pvals, method.str="Tukey HSD")
}


.newmankeuls <- function(center, means, n, MSE, df, conf.level) {
  se <- sqrt((MSE/2) * outer(1/n, 1/n, "+"))
  keep <- lower.tri(se)
  
  nmean <- (abs(outer(rank(means), rank(means), "-")) + 1)[keep]
  
  width <- qtukey(conf.level, nmean, df) * se[keep]
  est <- center / se[keep]
  pvals <- ptukey(abs(est), nmean, df, lower.tail = FALSE)
  
  list(width=width, est=est, pvals=pvals, method.str="Newman-Keuls")
}


.duncan <- function(center, means, n, MSE, df, conf.level) {
  se <- sqrt((MSE/2) * outer(1/n, 1/n, "+"))
  keep <- lower.tri(se)
  
  nmean <- (abs(outer(rank(means), rank(means), "-")) + 1)[keep]
  
  width <- qtukey(conf.level^(nmean-1), nmean, df) * se[keep]
  est <- center / se[keep]
  pvals <- 1 - (1 - ptukey(abs(est), nmean, df, lower.tail = FALSE))^(1/(nmean - 1))
  
  list(width=width, est=est, pvals=pvals,
       method.str="Duncan's new multiple range test")
}


.scheffe <- function(center, n, MSE, df, k, conf.level) {
  se <- sqrt(MSE * outer(1/n, 1/n, "+"))
  keep <- lower.tri(se)
  
  Fcrit <- qf(conf.level, k - 1, df)
  width <- sqrt((k - 1) * Fcrit) * se[keep]
  est <- center / se[keep]
  pvals <- pf(est^2/(k - 1), k - 1, df, lower.tail = FALSE)
  
  list(width=width, est=est, pvals=pvals, method.str="Scheff\u00e9")
}


