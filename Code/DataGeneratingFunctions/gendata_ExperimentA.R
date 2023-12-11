gendata_ExperimentA = function(n = 100, p = 20, m = 50, r = 0, depX = TRUE) {
  
  # Recall the generation of Y=Qin in example 6.1.1: Qin = Vmu + Vsigma%*%t(Ninvpgrid)
  mseq = seq(1 / (2 * m), 1 - 1 / (2 * m), 1 / m)
  
  Ninvpgrid = qnorm(mseq)
  x = 2 * pnorm(matrix(rnorm(n * p), n, p) %*% chol(r^(as.matrix(dist(1:p))))) - 1
  
  # Set initial parameters:
  Vmu0 = 0       # 0
  Vsigma0 = 3    # 3
  Vbeta = 3/4    # 3/4
  # Vbeta_new = 0.5
  Vgamma = 1     # 1
  Vnu1 = 1       # 1
  Vnu2 = 0.5     # 0.5
  
  ##  Vmu = Vmu0 + Vbeta * (x[ , 2] + x[ , 3]) + Vbeta_new * x[ , 4] + rnorm(n) * sqrt(Vnu1)
  if(depX){
    
    Vmu = Vmu0 + Vbeta * (x[ , 2] + x[ , 3]) + rnorm(n) * sqrt(Vnu1)
    GamP1 = (Vsigma0 + Vgamma * x[ , 1])^2 / Vnu2
    GamP2 = Vnu2 / (Vsigma0 + Vgamma * x[ , 1])
    
  } else {
    
    Vmu = Vmu0 + rnorm(n) * sqrt(Vnu1)
    GamP1 = (Vsigma0 + Vgamma)^2 / Vnu2
    GamP2 = Vnu2 / (Vsigma0 + Vgamma)
    
  }
  
  Vsigma = rgamma(n, shape = GamP1, scale = GamP2)
  # for (j in 1:n) Vsigma[j] = rgamma(1, shape = GamP1[j], scale = GamP2[j])
  
  Qin = tcrossprod(Vsigma, Ninvpgrid) + Vmu
  
  return(list('X' = x, 'Y' = Qin))
  
}
