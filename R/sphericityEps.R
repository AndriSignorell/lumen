
#' Greenhouse-Geisser And Huynh-Feldt Epsilons
#' 
#' Computes the Greenhouse-Geisser and Huynh-Feldt epsilon correction factors 
#' for assessing and correcting for violations of the sphericity assumption 
#' in repeated measures ANOVA.
#'  
#' @param S pxp covariance matrix
#' @param p dimension of observation vectors
#' @param g number of groups
#' @param n number of subjects
#' @param method a character string specifying which epsilon to return,
#'   must be one of \code{"both"} (default), \code{"GG"} for 
#'   Greenhouse-Geisser, or \code{"HF"} for Huynh-Feldt.
#'    
#' @return a numeric value

#' @note
#' Based on code by Hans Rudolf Roth.
#' 
#' @seealso \code{\link{aov}} 
#' 
#' @references Vonesh, E.F., Chinchilli, V.M. (1997) \emph{Linear and Nonlinear
#' Models for the Analysis of Repeated Measurements} Marcel Dekker, New York,
#' p.84-86
#' 
#' Crowder, M.J., Hand, D.J. (1990) \emph{Analysis of Repeated Measures}.
#' Chapman & Hall, London, p.54-55 
#' 
#' @examples
#' 
#' ## find!
#' 
#' 


#' @export
sphericityEps <- function(S, p, g, n, method = c("both", "GG", "HF")) {

  ## Purpose: calculates the Greenhouse-Geisser and Huynh-Feldt epsilons
  ## -------------------------------------------------------------------
  ## Arguments: S pxp covariance matrix
  ##            p dimension of observation vectors
  ##            g number of groups
  ##            n number of subjects
  
  ## Lit:    E.F. Vonesh + V.M. Chinchilli (1997), p.84-86
  ##         M.J. Crowder and D.J. Hand (1990), p.54-55
  
  ## Author: H.-R. Roth
  ## Date:   23.07.2002
  ## -------------------------------------------------------------------
  
  method <- match.arg(method)
  
  # U is a matrix of (p-1) orthonormal contrasts
  U <- t(cbind(diag(p - 1), 0) - outer(1:(p - 1), 1:p, "<") / ((p - 1):1))
  a <- 1 / sqrt(colSums(U^2))
  U <- U %*% diag(a)
  V <- t(U) %*% S %*% U
  
  e <- (sum(diag(V)))^2 / sum(diag(V %*% V)) / (p - 1)
  
  GG <- e
  HF <- min(1, (n * (p - 1) * e - 2) / ((p - 1) * (n - g - (p - 1) * e)))
  
  switch(method,
         "GG"   = c("GG-epsilon" = GG),
         "HF"   = c("HF-epsilon" = HF),
         "both" = c("GG-epsilon" = GG, "HF-epsilon" = HF)
  )
}


