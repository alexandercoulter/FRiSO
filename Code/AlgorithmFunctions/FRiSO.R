#' Fréchet Ridge Selection Operator (FRiSO)
#'
#' @param X An n by p matrix column-wise consisting of predictor vectors
#' @param Y An n by m matrix row-wise consisting of empirical quantile functions, each evaluated on an equispaced m-grid on (0, 1)
#' @param lower Lower box constraint on quantile functions; default is -Inf, i.e. unbounded
#' @param upper Upper box constraint on quantile functions; default is +Inf, i.e. unbounded
#' @param tauseq A p-length vector containing 'tau' values at which to solve FRiSO problem
#' @param lambda_init An optional p-length vector giving the initial 'allowance vector' lambda for FRiSO algorithm; will be scaled to sit on tau-simplex for first entry of tauseq
#' @param eps A non-negative error tolerance parameter
#' @param nudge A non-negative numeric value to offset warm starts to avoid spurious boundary values
#' @param alpha A non-negative dampening parameter (used for 'GSS' and 'SSCG' methods)
#' @param indicator_columns An optional vector of integers, from 1 to p, identifying indicator/binary columns in X
#' @param method String to determine method used to solve FRiSO problem, one of 'GSS' (geodesic short-step), 'FW' (Frank-Wolfe algorithm), or 'SSCG' (short-step conditional gradient)
#'
#' @return A p by length(tauseq) matrix column-wise containing fitted 'allowance vectors' lambda per 'tau' in tauseq.
#' @export
#'
#' @examples
FRiSO = function(X,
                 Y,
                 lower = -Inf,
                 upper = Inf,
                 tauseq,
                 lambda_init = NULL,
                 eps = 1e-5,
                 nudge = 0,
                 alpha = 0.9,
                 doNotScale = NULL,
                 method = c('GSS', 'FW', 'SSCG', 'SSAS'),
                 nesterov = FALSE,
                 maxIter = 1000,
                 verbose = FALSE){
  
  # Grab dimensions:
  n = nrow(X)
  m = ncol(Y)
  p = ncol(X)
  
  # Compatibility checks:
  if(n != nrow(Y)) stop("\'X\' and \'Y\' should have the same number of rows.")
  if(lower > upper) stop('Lower bound should be strictly less than upper bound.')
  
  if(is.null(lambda_init)) lambda_init = rep(tauseq[1] / p, p)
  if(length(lambda_init) != p) stop("\'lambda_init\' length and ncol(X) should match.")
  if(eps < 0) stop("\'eps\' should be non-negative.")
  
  # Center and scale inputs:
  # ALL NEEDED ROOT-n CALCULATIONS WILL BE DONE IN GSS FUNCTION
  Xc = X - rep(1, n) %*% crossprod(rep(1 / n, n), X)
  Xc_sd = sqrt(colMeans(Xc * Xc))
  Xc_sd[abs(Xc_sd) < 1e-10] = Inf
  # Xc = Xc / tcrossprod(rep(1, n), Xc_sd)
  if(TRUE){
    if(!is.null(doNotScale)){
      if((max(doNotScale) > p) | (min(doNotScale) < 1)) stop('doNotScale out of bounds of [1:p].')
      Xc[ , -doNotScale] = Xc[ , -doNotScale] / tcrossprod(rep(1, n), Xc_sd[-doNotScale])
    } else {
      Xc = Xc / tcrossprod(rep(1, n), Xc_sd)
    }
  }
  
  method = match.arg(method)
  if(method == 'GSS'){
    
    if(verbose){
      
      friso = FRiSO_GSS_Verbose(X = Xc,
                                Y = Y,
                                gamma_init = sqrt(lambda_init),
                                tauseq = tauseq,
                                lower = lower,
                                upper = upper,
                                alpha = alpha,
                                nudge = nudge * tauseq,
                                eps = eps,
                                nesterov = nesterov,
                                J = maxIter)
      
      all_LAMBDA = friso$All_LAMBDA
      LAMBDA = friso$LAMBDA
      objectives = friso$Objectives / (n * m)
      errors = friso$Errors
      return(list('LAMBDA' = LAMBDA, 'objectives' = objectives, 'errors' = errors, 'all_LAMBDA' = all_LAMBDA))
      
    } else {
      
      friso = FRiSO_GSS(X = Xc,
                        Y = Y,
                        gamma_init = sqrt(lambda_init),
                        tauseq = tauseq,
                        lower = lower,
                        upper = upper,
                        alpha = alpha,
                        nudge = nudge * tauseq,
                        eps = eps,
                        nesterov = nesterov,
                        J = maxIter)
      return(list('LAMBDA' = friso$LAMBDA))
      
    }
    
  } else if(method == 'FW'){
    
    friso = FRiSO_FW(X = Xc,
                     Y = Y,
                     lambda_init = lambda_init,
                     tauseq = tauseq,
                     lower = lower,
                     upper = upper,
                     eps = eps,
                     J = maxIter)
    all_LAMBDA = friso$All_LAMBDA
    LAMBDA = friso$LAMBDA
    objectives = friso$Objectives / (n * m)
    errors = friso$Errors
    return(list('LAMBDA' = LAMBDA, 'objectives' = objectives, 'errors' = errors, 'all_LAMBDA' = all_LAMBDA))
    
  } else if(method == 'SSCG'){
    
    friso = FRiSO_SSCG(X = Xc,
                       Y = Y,
                       lambda_init = lambda_init,
                       tauseq = tauseq,
                       lower = lower,
                       upper = upper,
                       a = alpha,
                       eps = eps,
                       J = maxIter)
    all_LAMBDA = friso$All_LAMBDA
    LAMBDA = friso$LAMBDA
    objectives = friso$Objectives / (n * m)
    errors = friso$Errors
    return(list('LAMBDA' = LAMBDA, 'objectives' = objectives, 'errors' = errors, 'all_LAMBDA' = all_LAMBDA))
    
  } else if(method == 'SSAS'){
    
    friso = FRiSO_SSAS(X = Xc,
                       Y = Y,
                       lambda_init = lambda_init,
                       tauseq = tauseq,
                       lower = lower,
                       upper = upper,
                       eps = eps,
                       J = maxIter)
    all_LAMBDA = friso$All_LAMBDA
    LAMBDA = friso$LAMBDA
    objectives = friso$Objectives / (n * m)
    errors = friso$Errors
    return(list('LAMBDA' = LAMBDA, 'objectives' = objectives, 'errors' = errors, 'all_LAMBDA' = all_LAMBDA))
    
  }

}
