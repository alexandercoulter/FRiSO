gendata_ExperimentB = function(n = 100, p = 20, m = 50, r = 0, depX = TRUE){
  
  U = chol(r^(as.matrix(dist(1:p))))
  # X = matrix(rnorm(n * p), n, p) %*% U
  X = 2 * pnorm(matrix(rnorm(n * p), n, p) %*% U) - 1
  
  # Control variables:
  # Baseline log(size):
  Vlogr0 = log(10)
  
  # Baseline logit(prob) of success:
  Vlogitp0 = log(0.5 / (1 - 0.5))
  
  # Baseline logit(prob) of zero:
  Vz0 = log(0.2 / (1 - 0.2))
  
  # log(size) adjustment due to X2, X3 covariates:
  Vbeta = 0.2
  
  # logit(prob) adjustment from X1 covariate:
  Vgamma = 0.2
  
  # logit(z) adjustment due to X4 covariate:
  Vpsi = 0.4
  
  # Baseline log(size) variability:
  Vnu1 = 0.15
  
  # Baseline logit(prob) variability:
  Vnu2 = 0.15
  
  # Baseline logit(z) variability:
  Vnu3 = 0.1
  
  if(depX){
    
    # SIZE:
    Vr = exp(rnorm(n, Vlogr0 + Vbeta * (X[ , 2] + X[ , 3]), Vnu1))
    
    # PROB SUCCESS:
    Vp = 1 / (1 + exp(-1 * rnorm(n, Vlogitp0 + Vgamma * X[ , 1], Vnu2)))
    
    # PROB ZERO:
    Vz = 1 / (1 + exp(-1 * rnorm(n, Vz0 + Vpsi * X[ , 4], Vnu3)))
    
  } else {
    
    # SIZE:
    Vr = exp(rnorm(n, Vlogr0, Vnu1))
    
    # PROB SUCCESS:
    Vp = 1 / (1 + exp(-1 * rnorm(n, Vlogitp0, Vnu2)))
    
    # PROB ZERO:
    Vz = 1 / (1 + exp(-1 * rnorm(n, Vz0, Vnu3)))
    
  }

  # Generate quantile functions:
  mseq = seq(1 / (2 * m), 1 - 1 / (2 * m), 1 / m)
  mseq = rep(mseq, each = n)
  mseq = pmax((mseq - rep(Vz, m)) / (1 - rep(Vz, m)), 0)
  Q = matrix(qnbinom(mseq, size = Vr, prob = Vp), n, m)
  
  return(list('X' = X, 'Y' = Q))
  # return(list('X' = X, 'Y' = Q, 'Xout' = Xout))
  
}
