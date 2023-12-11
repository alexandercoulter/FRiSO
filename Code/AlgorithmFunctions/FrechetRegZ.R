#' Title
#'
#' @param X 
#' @param Y 
#' @param Z 
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
FrechetRegZ = function(X,
                       Y,
                       Z = NULL,
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
  if(is.null(Z)) Z = X
  if(ncol(X) != ncol(Z)) stop('X and Z must have the same number of columns.')
  
  # Get dimensions:
  n = nrow(X)
  nz = nrow(Z)
  m = ncol(Y)
  p = ncol(X)
  
  # Center inputs:
  Xc = X - tcrossprod(rep(1, n), colMeans(X))
  Zc = Z - tcrossprod(rep(1, nz), colMeans(X))
  
  Xc_sd = sqrt(colMeans(Xc * Xc))
  Xc_sd[abs(Xc_sd) < 1e-10] = Inf
  # Xc = Xc / tcrossprod(rep(1, n), Xc_sd)
  if(TRUE){
    if(!is.null(doNotScale)){
      if((max(doNotScale) > p) | (min(doNotScale) < 1)) stop('doNotScale out of bounds of [1:p].')
      Zc[ , -doNotScale] = Zc[ , -doNotScale] / tcrossprod(rep(1, nz), Xc_sd[-doNotScale])
      Xc[ , -doNotScale] = Xc[ , -doNotScale] / tcrossprod(rep(1, n), Xc_sd[-doNotScale])
    } else {
      Zc = Zc / tcrossprod(rep(1, nz), Xc_sd)
      Xc = Xc / tcrossprod(rep(1, n), Xc_sd)
    }
  }
  
  Xc = Xc / sqrt(n)
  Zc = Zc / sqrt(n)
  
  # Calculate Yhat:
  # If 'lambda' is specified, calculate with lambda:
  if(!is.null(lambda)){
    
    # Dimension checking for method:
    if(p >= n){
      
      XDXt = tcrossprod(Xc * tcrossprod(rep(1, n), lambda), Xc)
      diag(XDXt) = diag(XDXt) + 1
      ZDXt = tcrossprod(Zc * tcrossprod(rep(1, nz), lambda), Xc)
      Yhat = rep(1, nz) %*% crossprod(rep(1 / n, n), Y) + ZDXt %*% solve(XDXt, Y)
      
    } else {
      
      XD = Xc * tcrossprod(rep(1, n), sqrt(lambda))
      ZD = Zc * tcrossprod(rep(1, nz), sqrt(lambda))
      G = crossprod(XD)
      diag(G) = diag(G) + 1
      Yhat = rep(1, nz) %*% crossprod(rep(1 / n, n), Y) + ZD %*% solve(G, crossprod(XD, Y))
      
    }
    
  } else { # If 'lambda' is not specified:
    
    invSigma = tryCatch( { solve(crossprod(Xc)) }, error = function(e){ MASS::ginv(crossprod(Xc)) } )
    Yhat = rep(1, nz) %*% crossprod(rep(1 / n, n), Y) + Zc %*% invSigma %*% crossprod(Xc, Y)
    
  }
  
  # Calculate monotone Q:
  Q = Qhat(Yhat, eps = eps)
  
  # Apply box constraints:
  Q[Q < lower] = lower
  Q[Q > upper] = upper
  
  # Return result:
  return(Q)
  
}
