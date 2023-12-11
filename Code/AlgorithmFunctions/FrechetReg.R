#' Title
#'
#' @param X 
#' @param Y 
#' @param lower 
#' @param upper 
#' @param lambda 
#' @param eps 
#' @param indicator_columns 
#'
#' @return
#' @export
#'
#' @examples
FrechetReg = function(X,
                      Y,
                      lower = -Inf,
                      upper = Inf,
                      lambda = NULL,
                      eps = 1e-5,
                      doNotScale = NULL){
  
  # Compatibility checks:
  if(nrow(X) != nrow(Y)) stop('Y and X should have the same number of rows.')
  if(!is.null(lambda)){
    if(ncol(X) != length(lambda)) stop('lambda should have the same length as X has rows.')
  }
  if(lower > upper) stop('Lower bound should be strictly less than upper bound.')
  
  # Get dimensions:
  n = nrow(X)
  m = ncol(Y)
  p = ncol(X)
  
  # Center and scale inputs:
  # ALL NEEDED ROOT-n CALCULATIONS WILL BE DONE IN GRADIENT DESCENT FUNCTIONS
  Xc = X - tcrossprod(rep(1, n), colMeans(X))
  Xc_sd = sqrt(colMeans(Xc * Xc))
  Xc_sd[abs(Xc_sd) < 1e-10] = Inf
  if(TRUE){
    if(!is.null(doNotScale)){
      if(max(doNotScale) > p) stop('Indicator columns out of range')
      Xc[ , -doNotScale] = Xc[ , -doNotScale] / tcrossprod(rep(1, n), Xc_sd[-doNotScale])
    } else {
      Xc = Xc / tcrossprod(rep(1, n), Xc_sd)
    }
  }
  Xc = Xc / sqrt(n)
  
  # Calculate unconstrained solution, Yhat:
  if(!is.null(lambda)){
    
    XcD = Xc * tcrossprod(rep(1, n), lambda)
    
    # If p >= n, then it's faster to do the "fast ridge" approach:
    if(p >= n){
      
      G = tcrossprod(XcD, Xc)
      diag(G) = diag(G) + 1
      Yhat = rep(1, n) %*% crossprod(rep(1 / n, n), Y) + Y - solve(G, Y)
      
    } else {
      
      G = crossprod(XcD, Xc)
      diag(G) = diag(G) + 1
      Yhat = rep(1, n) %*% crossprod(rep(1 / n, n), Y) + Xc %*% solve(G, crossprod(XcD, Y))
      
    }

  } else {
    
    U = svd(Xc)$u[ , -n]
    Yhat = rep(1, n) %*% crossprod(rep(1 / n, n), Y) + U %*% crossprod(U, Y)
    
  }
  
  # Calculate monotone Qhat:
  Q = Qhat(Yhat, eps)
  
  # Apply box constraints:
  Q[Q < lower] = lower
  Q[Q > upper] = upper
  
  # Return result:
  return(Q)
  
}
