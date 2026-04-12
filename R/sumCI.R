
#' Add Up Partial Confidence Intervals to a Total CI
#' 
#' Starting with a response variable that obtains different 
#' confidence intervals (CI) when calculated with different 
#' explanatory variables, all the values of the response variable 
#' should be added up. This function returns the CI for the sum.
#' 
#' @param x a matrix with 3 columns, containing the estimate in the first column 
#' followed by the lower and the upper confidence interval .

#' @return a vector with the sum and the lower, upper confidence 
#' bound of the confidence interval

#' @seealso \code{\link{binomCI}},

#' @references \url{https://stats.stackexchange.com/questions/223924/how-to-add-up-partial-confidence-intervals-to-create-a-total-confidence-interval}

#' @examples
#' x <- do.call(rbind, 
#'              tapply(bedrock::d.pizza$delivery_min, 
#'                     bedrock::d.pizza$area, meanCI))
#' sumCI(x)
#'


#' @export 
sumCI <- function(x){
  
  # https://stats.stackexchange.com/questions/223924/how-to-add-up-partial-confidence-intervals-to-create-a-total-confidence-interval
  # x: matrix with columns est, lci, uci
  
  # --- checks ---
  if (!is.matrix(x) || ncol(x) != 3) {
    stop("`x` must be a matrix with 3 columns: est, lci, uci.", call. = FALSE)
  }
  
  # --- extract ---
  est <- x[, 1]
  lci <- x[, 2]
  uci <- x[, 3]
  
  # --- half width ---
  hw <- (uci - lci) / 2
  
  # --- combine ---
  hw_sum <- sqrt(sum(hw^2))
  est_sum <- sum(est)
  
  c(
    sum = est_sum,
    lci = est_sum - hw_sum,
    uci = est_sum + hw_sum
  )
  
}

