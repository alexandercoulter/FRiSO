#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat Eta(const arma::mat& Yhat,
              const arma::mat& L,
              const arma::mat& U,
              const double& eps = 0.00001){
  
  // Initialize objects for gradient ascent:
  double n = Yhat.n_rows;
  double m = Yhat.n_cols;
  
  // Initialize Lagrange multiplier H:
  arma::mat H(n, m + 1, arma::fill::zeros);
  
  // Create constraint matrix:
  arma::mat C = arma::join_rows(L, Yhat) - arma::join_rows(Yhat, U);
  
  // Find rows where constraints are violated:
  arma::colvec maxC = max(C, 1);
  arma::uvec posC = find(maxC > 0);
  
  // If no rows have constraints violated, return H as-is:
  if(posC.n_elem == 0){
    
    return(H);
    
  } else {
    
    // Lagrange method:
    
    // Initialize - extract only those rows which will be worked on:
    arma::mat Hwork = H.rows(posC);
    arma::mat oldHwork = Hwork;
    arma::mat Cwork = C.rows(posC);
    Cwork *= 0.5;
    double error = eps + 1.0;
    
    // Go through loop:
    while(error > eps){
      
      // Set old H:
      oldHwork = Hwork;
      
      // Perform gradient ascent step
      Hwork.head_cols(m) = oldHwork.tail_cols(m);
      Hwork.tail_cols(m) += oldHwork.head_cols(m);
      Hwork.head_cols(1) += oldHwork.head_cols(1);
      Hwork *= 0.5;
      Hwork += Cwork;
      
      // Project to positive part:
      Hwork.transform( [](double val) { return (val <= 0.0) ? 0.0 : val; } );
      
      // Calculate error term:
      error = (abs(oldHwork - Hwork)).max();
      
    }
    
    H.rows(posC) = Hwork;
    return(H);
    
  }
  
}

// [[Rcpp::export]]
arma::mat add_I(const arma::mat& M,
                const double& a = 1){
  
  arma::mat Mnew = M;
  Mnew.diag() += a;
  return(Mnew);
  
}

// [[Rcpp::export]]
arma::mat XDXt(const arma::mat& X,
               const arma::colvec& d){
  
  arma::mat Xd = X;
  Xd.each_row() %= d.t();
  return(X * Xd.t());
  
}

// [[Rcpp::export]]
Rcpp::List FRiSO_GSS_Verbose(const arma::mat& X,
                             const arma::mat& Y,
                             const arma::colvec& gamma_init,
                             const arma::colvec& tauseq,
                             const arma::colvec& nudge,
                             const double& lower = 0,
                             const double& upper = 1,
                             const double& alpha = 1,
                             const double& eps = 0.0001,
                             const double& max_theta = 0.78539816,
                             const bool& nesterov = false,
                             const int& J = 1000){
  
  // Grab dimensions:
  double n = Y.n_rows;
  double m = Y.n_cols;
  double p = X.n_cols;
  
  // Initialize common double values:
  double tau;       // (current) tau value
  double tau_root;  // square root of tau value
  double tau_sq;    // square of tau value
  double ell = 0.0; // for Nesterov acceleration
  double eta = 0.0; // for Nesterov acceleration
  double der1;      // d / d\theta
  double der2;      // d2 / d\theta2
  double normv;     // norm of negative tangent gradient
  double theta;     // angle to rotate
  
  // Initialize common integers:
  int n_tau = tauseq.n_elem; // number of tau values to evaluate at
  int c;                     // number of active constraints
  
  // Initialize common vectors for algorithm:
  arma::uvec w;                        // vector for defining A0
  arma::colvec gamma_new = gamma_init; // gamma vector for looping
  arma::colvec gamma_old = gamma_init; // other gamma vector for looping
  arma::colvec phi_new = gamma_init;   // vector for Nesterov acceleration
  arma::colvec phi_old = gamma_init;   // other vector for Nesterov acceleration
  arma::colvec u(p);                   // normal vector for rotations
  arma::colvec v(p);                   // tangent vector for rotations
  arma::colvec sum_C(n);               // vector samples w/active constraints
  arma::colvec dn(p);                  // diagonal of N matrix
  arma::colvec grad(p);                // gradient vector
  arma::colvec grad_tan(p);            // tangent component of gradient vector
  arma::colvec L(n, arma::fill::value(lower)); // lower box constraint vector
  arma::colvec U(n, arma::fill::value(upper)); // upper box constraint vector
  
  // Set up matrices and vectors for function:
  arma::mat Xn = X / sqrt(n); // scaled X matrix (i.e. tilX)
  arma::mat Yt = Y.t();       // transposed Y (calculate once)
  arma::mat Xnt = Xn.t();     // transposed scaled X matrix (calculate once)
  arma::mat H(n, m + 1);      // Lagrange multiplier (capital Eta)
  arma::mat Yhat(n, m);       // unconstrained quantile (embedded) solution
  arma::mat Q(n, m);          // constrained quantile (embedded) solution
  arma::mat E(n, m);          // residuals matrix
  arma::mat tilE(n, m);       // adjusted residuals matrix
  arma::mat C(n, m - 1);      // active constraint matrix
  arma::mat V(m, n);          // matrix for second derivative
  arma::mat W(m, n);          // other matrix for second derivative
  
  // Set up objects to store outputs:
  arma::mat GAMMA(p, tauseq.n_elem); // matrix to store all FRiSO solutions (as
                                     // gammas)
  arma::colvec objs(J);              // vector to store all objective function
                                     // values for the last tau in tauseq
  arma::colvec errors(J);            // vector to store all errors for the last
                                     // tau in tauseq
  arma::mat all_LAMBDA(p, J);        // matrix to store all calculated FRiSO
                                     // solutions (as lambdas)

  // Two routes depending on how p and n compare:
  if(p > (1.2 * n)){
    
    // Initialize for (p > 1.2 * n) case:
    arma::mat YbarY = arma::mat(n, 1, arma::fill::ones) * (arma::mat(1, n, arma::fill::value(1 / n)) * Y) + Y;
    arma::mat Ginv(n, n);   // Xn D Xn' + I
    arma::mat GX(n, p);     // solve(Xn D Xn' + I, Xn)
    arma::mat XtG(p, n);    // solve(Xn D Xn' + I, Xn)'
    arma::mat GY(n, m);     // solve(Xn D Xn' + I, Y)
    arma::mat DuvXtG(p, n); // solve(Xn D Xn' + I, Xn Du Dv)'
    arma::mat YtGX(m, p);   // Y' solve(Xn D Xn' + I, Xn)
    
    // Loop through tau:
    for(int t = 0; t < n_tau; t++){
      
      /////////////////////////////////////////////////////////////////////////
      // Setup with new tau:
      //
      // Set current tau value:
      tau = tauseq(t);
      tau_root = sqrt(tau);
      tau_sq = tau * tau;
      
      // Update gamma value, 'nudging' it off the boundary:
      gamma_new += nudge(t);
      gamma_new *= tau_root / arma::norm(gamma_new);
      gamma_old = gamma_new;
      u = gamma_old / tau_root;
      
      if(nesterov){
        
        phi_new = gamma_new;
        phi_old = gamma_old;
        
      }
      
      // Reset ell and eta:
      ell = 0;
      eta = 0;
      
      // While loop:
      int j = 0;
      bool continueLoop = true;
      all_LAMBDA.zeros();
      
      while(continueLoop & (j < J)){
        
        ///////////////////////////////////////////////////////////////////////////
        // Set old values:
        gamma_old = gamma_new;
        phi_old = phi_new;
        u = gamma_old / tau_root;
        
        ///////////////////////////////////////////////////////////////////////////
        // Global solution:
        Ginv = add_I(XDXt(Xn, gamma_old % gamma_old));
        Yhat = YbarY - arma::solve(Ginv, Y, arma::solve_opts::likely_sympd);
        
        // Lagrange multiplier and constrained solution:
        H = Eta(Yhat, L, U, 0.00001);
        Q = Yhat + (H.head_cols(m) - H.tail_cols(m));
        Q.clamp(lower, upper);

        // Residuals matrix:
        E = Q - Y;
        
        // Objective function value:
        objs(j) = accu(E % E);
        
        // Constraint matrix:
        C = H;
        C.transform( [](double val) { return (val > 0) ? 1.0 : 0.0; } );
        
        ////////////////////////////////////////////////////////////////////////
        // Step 1: Calculate GX:
        GX = arma::solve(Ginv, Xn, arma::solve_opts::likely_sympd);
        XtG = GX.t();
        
        // Step 2: Calculate YtGX:
        YtGX = Y.t() * GX;
        
        // Step 3: Loop through i = {1 ... n} to calculate tilE:
        tilE = E;
        sum_C = sum(C, 1);
        if(sum(sum_C) == 0){}else{
          for(int i = 0; i < n; i++){
            
            // If sum_C(i) is non-zero, we need to adjust this row of Eadj:
            if(sum_C(i) > 0){

              // Save sum_C(i):
              c = sum_C(i);
              
              // Get indices where C(i, ) is non-zero:
              w = find(C.row(i)) + arma::linspace<arma::uvec>(0, c - 1, c) * m;
              
              // Create A0 object:
              arma::mat A0(m, c, arma::fill::zeros);
              A0(w(arma::find(w < (m * c)))) += 1;
              A0(w(arma::find(w > 0)) - 1) += -1;
              
              // Adjust this row:
              tilE.row(i) -= (tilE.row(i) * A0) * arma::solve(A0.t() * A0, A0.t(), arma::solve_opts::likely_sympd);
              
            }
            
          }
        }
        
        // Step 4: Calculate dn:
        dn = sum(GX % (tilE * YtGX), 0).t();
        
        // Step 4.25: Calculate gradient:
        grad = 4 * gamma_old % dn;
        
        // Step 4.5: Calculate tangential gradient:
        grad_tan = grad - u * arma::dot(u, grad);
        
        // Step 4.75: Calculate v:
        normv = arma::norm(grad_tan);
        v = -grad_tan / normv;
        
        // Step 5: Calculate GXDuDv.t():
        DuvXtG = XtG;
        DuvXtG.each_col() %= (u % v);
        
        // Step 6: Calculate V:
        V = YtGX * DuvXtG;
        
        // Step 7: Calculate W:
        W = V * (Xn * DuvXtG);
        
        // Step 8: Calculate sum(W % tilE):
        der2 = -16 * tau_sq * arma::accu(W % tilE.t());
        
        // Step 9: Loop through i = {1 ... n} to calculate tilV and sum(V % tilV):
        if(sum(sum_C) == 0){}else{
          for(int i = 0; i < n; i++){
            
            // If sum_C(i) is non-zero, we need to adjust this row of Eadj:
            if(sum_C(i) > 0){
              
              // Save sum_C(i):
              c = sum_C(i);
              
              // Get indices where C(i, ) is non-zero:
              w = find(C.row(i)) + arma::linspace<arma::uvec>(0, c - 1, c) * m;
              
              // Create A0 object:
              arma::mat A0(m, c, arma::fill::zeros);
              A0(w(arma::find(w < (m * c)))) += 1;
              A0(w(arma::find(w > 0)) - 1) += -1;
              
              // Adjust this row:
              V.col(i) -= A0 * arma::solve(A0.t() * A0, A0.t() * V.col(i), arma::solve_opts::likely_sympd);
              
            }
            
          }
        }
        
        der2 += 8 * tau_sq * arma::accu(V % V);
        
        //////////////////////////////////////////////////////////////////////
        // Finish first derivative:
        der1 = -tau_root * normv;
        
        // Finish second derivative:
        der2 += 4 * tau * sum(dn % (v % v - u % u));
        
        // Calculate angle theta:
        if(der2 == 0){
          
          theta = max_theta;
          
        } else {
          
          theta = alpha * abs(der1 / der2);
          if(theta > max_theta){
            
            theta = max_theta;
            
          }
          
        }
        
        // Optional Nesterov acceleration:
        if(nesterov){
          
          // Update the scaling factor:
          eta = 1 - ell;
          ell = 0.5 * (1 + sqrt(1 + 4 * ell * ell));
          eta = eta / ell;
          
          // In the case we're in the first iteration, eta must be zero:
          if(j == 0) eta = 0;
          
          // Move phi:
          phi_new = abs(cos(theta) * gamma_old + sin(theta) * v * tau_root);
          phi_new = phi_new * tau_root / norm(phi_new);
          
          // Calculate new angle between phi_old and phi_new, representing magnitude of this step:
          theta = dot(phi_new, phi_old) / tau;
          if(theta > 1) theta = 1;
          theta = acos(theta);
          
          // Calculate new "start" and "end" locations:
          u = phi_old / tau_root;
          v = phi_new - u * arma::dot(u, phi_new);
          v = v / norm(v);
          
          // Calculate gamma_new as the (1 - eta) accelerated shift from phi_old toward phi_new:
          gamma_new = tau_root * abs(cos(theta * (1 - eta)) * u + sin(theta * (1 - eta)) * v);
          gamma_new = gamma_new * tau_root / norm(gamma_new);
          
        } else {
          
          // Rotate gamma_old:
          gamma_new = abs(cos(theta) * gamma_old + sin(theta) * v * tau_root);
          gamma_new = gamma_new * tau_root / norm(gamma_new);
          
        }
        
        errors(j) = abs(gamma_new - gamma_old).max();
        
        // Check if gamma_new is close to gamma_old:
        if(errors(j) < eps){
          
          continueLoop = false;
          
        }
        
        // Increment counter:
        all_LAMBDA.col(j) = gamma_new % gamma_new;
        j++;
        
      }
      
      // Set the solution for the current tau:
      GAMMA.col(t) = gamma_new;
      
    }
    
  } else {
    
    // Initialize for (p <= 1.2 * n) case:
    arma::mat Ybar = arma::mat(n, 1, arma::fill::ones) * (arma::mat(1, n, arma::fill::value(1 / n)) * Y);
    arma::mat XtY = Xnt * Y;    // Xn' Y
    arma::mat Sigma = Xnt * Xn; // Xn' Xn (estimated cov(X))
    arma::mat Ftinv(p, p);      // Sigma D + I
    arma::mat FtXt(p, n);       // solve(Sigma D + I, Xn')
    arma::mat FtXtY(p, m);      // solve(Sigma D + I, Xn' Y)
    arma::mat DFtXtY(p, m);     // D solve(Sigma D + I, Xn' Y)
    arma::mat DuvFtXt(p, n);    // Du Dv solve(Sigma D + I, Xn')
    
    // Loop through tau:
    for(int t = 0; t < n_tau; t++){
      
      /////////////////////////////////////////////////////////////////////////
      // Setup with new tau:
      //
      // Set current tau value:
      tau = tauseq(t);
      tau_root = sqrt(tau);
      tau_sq = tau * tau;
      
      // Update gamma value, 'nudging' it off the boundary:
      gamma_new += nudge(t);
      gamma_new *= tau_root / arma::norm(gamma_new);
      gamma_old = gamma_new;
      u = gamma_old / tau_root;
      
      if(nesterov){
        
        phi_new = gamma_new;
        phi_old = gamma_old;
        
      }
      
      // Reset ell and eta:
      ell = 0.0;
      eta = 0.0;
      
      // While loop:
      int j = 0;
      bool continueLoop = true;
      all_LAMBDA.zeros();
      
      while(continueLoop & (j < J)){
        
        ///////////////////////////////////////////////////////////////////////////
        // Set old values:
        gamma_old = gamma_new;
        phi_old = phi_new;
        u = gamma_old / tau_root;
        
        ////////////////////////////////////////////////////////////////////////
        // Step 1: Calculate FtXt:
        Ftinv = Sigma;
        Ftinv.each_row() %= (gamma_old % gamma_old).t();
        Ftinv.diag() += 1;
        
        FtXt = arma::solve(Ftinv, Xnt, arma::solve_opts::fast);
        
        // Step 2: Calculate FtXtY:
        FtXtY = arma::solve(Ftinv, XtY, arma::solve_opts::fast);
        
        // Global solution:
        DFtXtY = FtXtY;
        DFtXtY.each_col() %= (gamma_old % gamma_old);
        Yhat = Ybar + Xn * DFtXtY;
        
        // Lagrange multiplier and constrained solution:
        H = Eta(Yhat, L, U, 0.00001);
        Q = Yhat + (H.head_cols(m) - H.tail_cols(m));
        Q.clamp(lower, upper);
        
        // Residuals matrix:
        E = Q - Y;
        
        // Objective function value:
        objs(j) = accu(E % E);
        
        // Constraint matrix:
        C = H;
        C.transform( [](double val) { return (val > 0) ? 1.0 : 0.0; } );
        
        // Step 3: Loop through i = {1 ... n} to calculate tilE:
        tilE = E;
        sum_C = sum(C, 1);
        if(sum(sum_C) == 0){}else{
          for(int i = 0; i < n; i++){
            
            // If sum_C(i) is non-zero, we need to adjust this row of Eadj:
            if(sum_C(i) > 0){
              
              // Save sum_C(i):
              c = sum_C(i);
              
              // Get indices where C(i, ) is non-zero:
              w = find(C.row(i)) + arma::linspace<arma::uvec>(0, c - 1, c) * m;
              
              // Create A0 object:
              arma::mat A0(m, c, arma::fill::zeros);
              A0(w(arma::find(w < (m * c)))) += 1;
              A0(w(arma::find(w > 0)) - 1) += -1;
              
              // Adjust this row:
              tilE.row(i) -= (tilE.row(i) * A0) * arma::solve(A0.t() * A0, A0.t(), arma::solve_opts::likely_sympd);
              
            }
            
          }
        }
        
        // Step 4: Calculate dn:
        dn = sum(FtXt % (FtXtY * tilE.t()), 1);
        
        // Step 4.25: Calculate gradient:
        grad = 4 * gamma_old % dn;
        
        // Step 4.5: Calculate tangential gradient:
        grad_tan = grad - u * arma::dot(u, grad);
        
        // Step 4.75: Calculate v:
        normv = arma::norm(grad_tan);
        v = -grad_tan / normv;
        
        // Step 5: Calculate DuvFtXt:
        DuvFtXt = FtXt;
        DuvFtXt.each_col() %= (u % v);
        
        // Step 6: Calculate V:
        V = FtXtY.t() * DuvFtXt;
        
        // Step 7: Calculate W:
        W = (V * Xn) * DuvFtXt;
        
        // Step 8: Calculate sum(W % tilE):
        der2 = -16 * tau_sq * arma::accu(W % tilE.t());
        
        // Step 9: Loop through i = {1 ... n} to calculate tilV and sum(V % tilV):
        if(sum(sum_C) == 0){}else{
          for(int i = 0; i < n; i++){
            
            // If sum_C(i) is non-zero, we need to adjust this row of Eadj:
            if(sum_C(i) > 0){
              
              // Save sum_C(i):
              c = sum_C(i);
              
              // Get indices where C(i, ) is non-zero:
              w = find(C.row(i)) + arma::linspace<arma::uvec>(0, c - 1, c) * m;
              
              // Create A0 object:
              arma::mat A0(m, c, arma::fill::zeros);
              A0(w(arma::find(w < (m * c)))) += 1;
              A0(w(arma::find(w > 0)) - 1) += -1;
              
              // Adjust this row:
              V.col(i) -= A0 * arma::solve(A0.t() * A0, A0.t() * V.col(i), arma::solve_opts::likely_sympd);
              
            }
            
          }
        }
        
        der2 += 8 * tau_sq * arma::accu(V % V);
        
        //////////////////////////////////////////////////////////////////////
        // Finish first derivative:
        der1 = -tau_root * normv;
        
        // Finish second derivative:
        der2 += 4 * tau * sum(dn % (v % v - u % u));
        
        // Calculate angle theta:
        if(der2 == 0){
          
          theta = max_theta;
          
        } else {
          
          theta = alpha * abs(der1 / der2);
          if(theta > max_theta){
            
            theta = max_theta;
            
          }
          
        }
        
        // Optional Nesterov acceleration:
        if(nesterov){
          
          // Update the scaling factor:
          eta = 1 - ell;
          ell = 0.5 * (1 + sqrt(1 + 4 * ell * ell));
          eta = eta / ell;
          
          // In the case we're in the first iteration, eta must be zero:
          if(j == 0) eta = 0;
          
          // Move phi:
          phi_new = abs(cos(theta) * gamma_old + sin(theta) * v * tau_root);
          phi_new = phi_new * tau_root / norm(phi_new);
          
          // Calculate new angle between phi_old and phi_new, representing magnitude of this step:
          theta = dot(phi_new, phi_old) / tau;
          if(theta > 1) theta = 1;
          theta = acos(theta);
          
          // Calculate new "start" and "end" locations:
          u = phi_old / tau_root;
          v = phi_new - u * arma::dot(u, phi_new);
          v = v / norm(v);
          
          // Calculate gamma_new as the (1 - eta) accelerated shift from phi_old toward phi_new:
          gamma_new = tau_root * abs(cos(theta * (1 - eta)) * u + sin(theta * (1 - eta)) * v);
          gamma_new = gamma_new * tau_root / norm(gamma_new);
          
        } else {
          
          // Rotate gamma_old:
          gamma_new = abs(cos(theta) * gamma_old + sin(theta) * v * tau_root);
          gamma_new = gamma_new * tau_root / norm(gamma_new);
          
        }
        
        errors(j) = abs(gamma_new - gamma_old).max();
        
        // Check if gamma_new is close to gamma_old:
        if(errors(j) < eps){
          
          continueLoop = false;
          
        }
        
        // Increment counter:
        all_LAMBDA.col(j) = gamma_new % gamma_new;
        j++;
        
      }
      
      // Set the solution for the current tau:
      GAMMA.col(t) = gamma_new;
      
    }
    
  }
  
  // Return output list:
  return Rcpp::List::create(Rcpp::Named("LAMBDA") = GAMMA % GAMMA,
                            Rcpp::Named("Objectives") = objs,
                            Rcpp::Named("Errors") = errors,
                            Rcpp::Named("All_LAMBDA") = all_LAMBDA);
  
}

// [[Rcpp::export]]
Rcpp::List FRiSO_GSS(const arma::mat& X,
                     const arma::mat& Y,
                     const arma::colvec& gamma_init,
                     const arma::colvec& tauseq,
                     const arma::colvec& nudge,
                     const double& lower = 0,
                     const double& upper = 1,
                     const double& alpha = 1,
                     const double& eps = 0.0001,
                     const double& max_theta = 0.78539816,
                     const bool& nesterov = false,
                     const int& J = 1000){
  
  // Grab dimensions:
  double n = Y.n_rows;
  double m = Y.n_cols;
  double p = X.n_cols;
  
  // Initialize common double values:
  double tau;       // (current) tau value
  double tau_root;  // square root of tau value
  double tau_sq;    // square of tau value
  double ell = 0.0; // for Nesterov acceleration
  double eta = 0.0; // for Nesterov acceleration
  double der1;      // d / d\theta
  double der2;      // d2 / d\theta2
  double normv;     // norm of negative tangent gradient
  double theta;     // angle to rotate
  
  // Initialize common integers:
  int n_tau = tauseq.n_elem; // number of tau values to evaluate at
  int c;                     // number of active constraints
  
  // Initialize common vectors for algorithm:
  arma::uvec w;                        // vector for defining A0
  arma::colvec gamma_new = gamma_init; // gamma vector for looping
  arma::colvec gamma_old = gamma_init; // other gamma vector for looping
  arma::colvec phi_new = gamma_init;   // vector for Nesterov acceleration
  arma::colvec phi_old = gamma_init;   // other vector for Nesterov acceleration
  arma::colvec u(p);                   // normal vector for rotations
  arma::colvec v(p);                   // tangent vector for rotations
  arma::colvec sum_C(n);               // vector samples w/active constraints
  arma::colvec dn(p);                  // diagonal of N matrix
  arma::colvec grad(p);                // gradient vector
  arma::colvec grad_tan(p);            // tangent component of gradient vector
  arma::colvec L(n, arma::fill::value(lower)); // lower box constraint vector
  arma::colvec U(n, arma::fill::value(upper)); // upper box constraint vector
  
  // Set up matrices and vectors for function:
  arma::mat Xn = X / sqrt(n); // scaled X matrix (i.e. tilX)
  arma::mat Yt = Y.t();       // transposed Y (calculate once)
  arma::mat Xnt = Xn.t();     // transposed scaled X matrix (calculate once)
  arma::mat H(n, m + 1);      // Lagrange multiplier (capital Eta)
  arma::mat Yhat(n, m);       // unconstrained quantile (embedded) solution
  arma::mat Q(n, m);          // constrained quantile (embedded) solution
  arma::mat E(n, m);          // residuals matrix
  arma::mat tilE(n, m);       // adjusted residuals matrix
  arma::mat C(n, m - 1);      // active constraint matrix
  arma::mat V(m, n);          // matrix for second derivative
  arma::mat W(m, n);          // other matrix for second derivative
  
  // Set up objects to store outputs:
  arma::mat GAMMA(p, tauseq.n_elem); // matrix to store all FRiSO solutions (as
                                     // gammas)
  
  // Two routes depending on how p and n compare:
  if(p > (1.2 * n)){
    
    // Initialize for (p > 1.2 * n) case:
    arma::mat YbarY = arma::mat(n, 1, arma::fill::ones) * (arma::mat(1, n, arma::fill::value(1 / n)) * Y) + Y;
    arma::mat Ginv(n, n);   // Xn D Xn' + I
    arma::mat GX(n, p);     // solve(Xn D Xn' + I, Xn)
    arma::mat XtG(p, n);    // solve(Xn D Xn' + I, Xn)'
    arma::mat GY(n, m);     // solve(Xn D Xn' + I, Y)
    arma::mat DuvXtG(p, n); // solve(Xn D Xn' + I, Xn Du Dv)'
    arma::mat YtGX(m, p);   // Y' solve(Xn D Xn' + I, Xn)
    
    // Loop through tau:
    for(int t = 0; t < n_tau; t++){
      
      /////////////////////////////////////////////////////////////////////////
      // Setup with new tau:
      //
      // Set current tau value:
      tau = tauseq(t);
      tau_root = sqrt(tau);
      tau_sq = tau * tau;
      
      // Update gamma value, 'nudging' it off the boundary:
      gamma_new += nudge(t);
      gamma_new *= tau_root / arma::norm(gamma_new);
      gamma_old = gamma_new;
      u = gamma_old / tau_root;
      
      if(nesterov){
        
        phi_new = gamma_new;
        phi_old = gamma_old;
        
      }
      
      // Reset ell and eta:
      ell = 0;
      eta = 0;
      
      // While loop:
      int j = 0;
      bool continueLoop = true;
      
      while(continueLoop & (j < J)){
        
        ///////////////////////////////////////////////////////////////////////////
        // Set old values:
        gamma_old = gamma_new;
        phi_old = phi_new;
        u = gamma_old / tau_root;
        
        ///////////////////////////////////////////////////////////////////////////
        // Global solution:
        Ginv = add_I(XDXt(Xn, gamma_old % gamma_old));
        Yhat = YbarY - arma::solve(Ginv, Y, arma::solve_opts::likely_sympd);
        
        // Lagrange multiplier and constrained solution:
        H = Eta(Yhat, L, U, 0.00001);
        Q = Yhat + (H.head_cols(m) - H.tail_cols(m));
        Q.clamp(lower, upper);
        
        // Residuals matrix:
        E = Q - Y;
        
        // Constraint matrix:
        C = H;
        C.transform( [](double val) { return (val > 0) ? 1.0 : 0.0; } );
        
        ////////////////////////////////////////////////////////////////////////
        // Step 1: Calculate GX:
        GX = arma::solve(Ginv, Xn, arma::solve_opts::likely_sympd);
        XtG = GX.t();
        
        // Step 2: Calculate YtGX:
        YtGX = Y.t() * GX;
        
        // Step 3: Loop through i = {1 ... n} to calculate tilE:
        tilE = E;
        sum_C = sum(C, 1);
        if(sum(sum_C) == 0){}else{
          for(int i = 0; i < n; i++){
            
            // If sum_C(i) is non-zero, we need to adjust this row of Eadj:
            if(sum_C(i) > 0){
              
              // Save sum_C(i):
              c = sum_C(i);
              
              // Get indices where C(i, ) is non-zero:
              w = find(C.row(i)) + arma::linspace<arma::uvec>(0, c - 1, c) * m;
              
              // Create A0 object:
              arma::mat A0(m, c, arma::fill::zeros);
              A0(w(arma::find(w < (m * c)))) += 1;
              A0(w(arma::find(w > 0)) - 1) += -1;
              
              // Adjust this row:
              tilE.row(i) -= (tilE.row(i) * A0) * arma::solve(A0.t() * A0, A0.t(), arma::solve_opts::likely_sympd);
              
            }
            
          }
        }
        
        // Step 4: Calculate dn:
        dn = sum(GX % (tilE * YtGX), 0).t();
        
        // Step 4.25: Calculate gradient:
        grad = 4 * gamma_old % dn;
        
        // Step 4.5: Calculate tangential gradient:
        grad_tan = grad - u * arma::dot(u, grad);
        
        // Step 4.75: Calculate v:
        normv = arma::norm(grad_tan);
        v = -grad_tan / normv;
        
        // Step 5: Calculate GXDuDv.t():
        DuvXtG = XtG;
        DuvXtG.each_col() %= (u % v);
        
        // Step 6: Calculate V:
        V = YtGX * DuvXtG;
        
        // Step 7: Calculate W:
        W = V * (Xn * DuvXtG);
        
        // Step 8: Calculate sum(W % tilE):
        der2 = -16 * tau_sq * arma::accu(W % tilE.t());
        
        // Step 9: Loop through i = {1 ... n} to calculate tilV and sum(V % tilV):
        if(sum(sum_C) == 0){}else{
          for(int i = 0; i < n; i++){
            
            // If sum_C(i) is non-zero, we need to adjust this row of Eadj:
            if(sum_C(i) > 0){
              
              // Save sum_C(i):
              c = sum_C(i);
              
              // Get indices where C(i, ) is non-zero:
              w = find(C.row(i)) + arma::linspace<arma::uvec>(0, c - 1, c) * m;
              
              // Create A0 object:
              arma::mat A0(m, c, arma::fill::zeros);
              A0(w(arma::find(w < (m * c)))) += 1;
              A0(w(arma::find(w > 0)) - 1) += -1;
              
              // Adjust this row:
              V.col(i) -= A0 * arma::solve(A0.t() * A0, A0.t() * V.col(i), arma::solve_opts::likely_sympd);
              
            }
            
          }
        }
        
        der2 += 8 * tau_sq * arma::accu(V % V);
        
        //////////////////////////////////////////////////////////////////////
        // Finish first derivative:
        der1 = -tau_root * normv;
        
        // Finish second derivative:
        der2 += 4 * tau * sum(dn % (v % v - u % u));
        
        // Calculate angle theta:
        if(der2 == 0){
          
          theta = max_theta;
          
        } else {
          
          theta = alpha * abs(der1 / der2);
          if(theta > max_theta){
            
            theta = max_theta;
            
          }
          
        }
        
        // Optional Nesterov acceleration:
        if(nesterov){
          
          // Update the scaling factor:
          eta = 1 - ell;
          ell = 0.5 * (1 + sqrt(1 + 4 * ell * ell));
          eta = eta / ell;
          
          // In the case we're in the first iteration, eta must be zero:
          if(j == 0) eta = 0;
          
          // Move phi:
          phi_new = abs(cos(theta) * gamma_old + sin(theta) * v * tau_root);
          phi_new = phi_new * tau_root / norm(phi_new);
          
          // Calculate new angle between phi_old and phi_new, representing magnitude of this step:
          theta = dot(phi_new, phi_old) / tau;
          if(theta > 1) theta = 1;
          theta = acos(theta);
          
          // Calculate new "start" and "end" locations:
          u = phi_old / tau_root;
          v = phi_new - u * arma::dot(u, phi_new);
          v = v / norm(v);
          
          // Calculate gamma_new as the (1 - eta) accelerated shift from phi_old toward phi_new:
          gamma_new = tau_root * abs(cos(theta * (1 - eta)) * u + sin(theta * (1 - eta)) * v);
          gamma_new = gamma_new * tau_root / norm(gamma_new);
          
        } else {
          
          // Rotate gamma_old:
          gamma_new = abs(cos(theta) * gamma_old + sin(theta) * v * tau_root);
          gamma_new = gamma_new * tau_root / norm(gamma_new);
          
        }
        
        // Check if gamma_new is close to gamma_old:
        if(abs(gamma_new - gamma_old).max() < eps){
          
          continueLoop = false;
          
        }
        
        // Increment counter:
        j++;
        
      }
      
      // Set the solution for the current tau:
      GAMMA.col(t) = gamma_new;
      
    }
    
  } else {
    
    // Initialize for (p <= 1.2 * n) case:
    arma::mat Ybar = arma::mat(n, 1, arma::fill::ones) * (arma::mat(1, n, arma::fill::value(1 / n)) * Y);
    arma::mat XtY = Xnt * Y;    // Xn' Y
    arma::mat Sigma = Xnt * Xn; // Xn' Xn (estimated cov(X))
    arma::mat Ftinv(p, p);      // Sigma D + I
    arma::mat FtXt(p, n);       // solve(Sigma D + I, Xn')
    arma::mat FtXtY(p, m);      // solve(Sigma D + I, Xn' Y)
    arma::mat DFtXtY(p, m);     // D solve(Sigma D + I, Xn' Y)
    arma::mat DuvFtXt(p, n);    // Du Dv solve(Sigma D + I, Xn')
    
    // Loop through tau:
    for(int t = 0; t < n_tau; t++){
      
      /////////////////////////////////////////////////////////////////////////
      // Setup with new tau:
      //
      // Set current tau value:
      tau = tauseq(t);
      tau_root = sqrt(tau);
      tau_sq = tau * tau;
      
      // Update gamma value, 'nudging' it off the boundary:
      gamma_new += nudge(t);
      gamma_new *= tau_root / arma::norm(gamma_new);
      gamma_old = gamma_new;
      u = gamma_old / tau_root;
      
      if(nesterov){
        
        phi_new = gamma_new;
        phi_old = gamma_old;
        
      }
      
      // Reset ell and eta:
      ell = 0.0;
      eta = 0.0;
      
      // While loop:
      int j = 0;
      bool continueLoop = true;
      
      while(continueLoop & (j < J)){
        
        ///////////////////////////////////////////////////////////////////////////
        // Set old values:
        gamma_old = gamma_new;
        phi_old = phi_new;
        u = gamma_old / tau_root;
        
        ////////////////////////////////////////////////////////////////////////
        // Step 1: Calculate FtXt:
        Ftinv = Sigma;
        Ftinv.each_row() %= (gamma_old % gamma_old).t();
        Ftinv.diag() += 1;
        
        FtXt = arma::solve(Ftinv, Xnt, arma::solve_opts::fast);
        
        // Step 2: Calculate FtXtY:
        FtXtY = arma::solve(Ftinv, XtY, arma::solve_opts::fast);
        
        // Global solution:
        DFtXtY = FtXtY;
        DFtXtY.each_col() %= (gamma_old % gamma_old);
        Yhat = Ybar + Xn * DFtXtY;
        
        // Lagrange multiplier and constrained solution:
        H = Eta(Yhat, L, U, 0.00001);
        Q = Yhat + (H.head_cols(m) - H.tail_cols(m));
        Q.clamp(lower, upper);
        
        // Residuals matrix:
        E = Q - Y;
        
        // Constraint matrix:
        C = H;
        C.transform( [](double val) { return (val > 0) ? 1.0 : 0.0; } );
        
        // Step 3: Loop through i = {1 ... n} to calculate tilE:
        tilE = E;
        sum_C = sum(C, 1);
        if(sum(sum_C) == 0){}else{
          for(int i = 0; i < n; i++){
            
            // If sum_C(i) is non-zero, we need to adjust this row of Eadj:
            if(sum_C(i) > 0){
              
              // Save sum_C(i):
              c = sum_C(i);
              
              // Get indices where C(i, ) is non-zero:
              w = find(C.row(i)) + arma::linspace<arma::uvec>(0, c - 1, c) * m;
              
              // Create A0 object:
              arma::mat A0(m, c, arma::fill::zeros);
              A0(w(arma::find(w < (m * c)))) += 1;
              A0(w(arma::find(w > 0)) - 1) += -1;
              
              // Adjust this row:
              tilE.row(i) -= (tilE.row(i) * A0) * arma::solve(A0.t() * A0, A0.t(), arma::solve_opts::likely_sympd);
              
            }
            
          }
        }
        
        // Step 4: Calculate dn:
        dn = sum(FtXt % (FtXtY * tilE.t()), 1);
        
        // Step 4.25: Calculate gradient:
        grad = 4 * gamma_old % dn;
        
        // Step 4.5: Calculate tangential gradient:
        grad_tan = grad - u * arma::dot(u, grad);
        
        // Step 4.75: Calculate v:
        normv = arma::norm(grad_tan);
        v = -grad_tan / normv;
        
        // Step 5: Calculate DuvFtXt:
        DuvFtXt = FtXt;
        DuvFtXt.each_col() %= (u % v);
        
        // Step 6: Calculate V:
        V = FtXtY.t() * DuvFtXt;
        
        // Step 7: Calculate W:
        W = (V * Xn) * DuvFtXt;
        
        // Step 8: Calculate sum(W % tilE):
        der2 = -16 * tau_sq * arma::accu(W % tilE.t());
        
        // Step 9: Loop through i = {1 ... n} to calculate tilV and sum(V % tilV):
        if(sum(sum_C) == 0){}else{
          for(int i = 0; i < n; i++){
            
            // If sum_C(i) is non-zero, we need to adjust this row of Eadj:
            if(sum_C(i) > 0){
              
              // Save sum_C(i):
              c = sum_C(i);
              
              // Get indices where C(i, ) is non-zero:
              w = find(C.row(i)) + arma::linspace<arma::uvec>(0, c - 1, c) * m;
              
              // Create A0 object:
              arma::mat A0(m, c, arma::fill::zeros);
              A0(w(arma::find(w < (m * c)))) += 1;
              A0(w(arma::find(w > 0)) - 1) += -1;
              
              // Adjust this row:
              V.col(i) -= A0 * arma::solve(A0.t() * A0, A0.t() * V.col(i), arma::solve_opts::likely_sympd);
              
            }
            
          }
        }
        
        der2 += 8 * tau_sq * arma::accu(V % V);
        
        //////////////////////////////////////////////////////////////////////
        // Finish first derivative:
        der1 = -tau_root * normv;
        
        // Finish second derivative:
        der2 += 4 * tau * sum(dn % (v % v - u % u));
        
        // Calculate angle theta:
        if(der2 == 0){
          
          theta = max_theta;
          
        } else {
          
          theta = alpha * abs(der1 / der2);
          if(theta > max_theta){
            
            theta = max_theta;
            
          }
          
        }
        
        // Optional Nesterov acceleration:
        if(nesterov){
          
          // Update the scaling factor:
          eta = 1 - ell;
          ell = 0.5 * (1 + sqrt(1 + 4 * ell * ell));
          eta = eta / ell;
          
          // In the case we're in the first iteration, eta must be zero:
          if(j == 0) eta = 0;
          
          // Move phi:
          phi_new = abs(cos(theta) * gamma_old + sin(theta) * v * tau_root);
          phi_new = phi_new * tau_root / norm(phi_new);
          
          // Calculate new angle between phi_old and phi_new, representing magnitude of this step:
          theta = dot(phi_new, phi_old) / tau;
          if(theta > 1) theta = 1;
          theta = acos(theta);
          
          // Calculate new "start" and "end" locations:
          u = phi_old / tau_root;
          v = phi_new - u * arma::dot(u, phi_new);
          v = v / norm(v);
          
          // Calculate gamma_new as the (1 - eta) accelerated shift from phi_old toward phi_new:
          gamma_new = tau_root * abs(cos(theta * (1 - eta)) * u + sin(theta * (1 - eta)) * v);
          gamma_new = gamma_new * tau_root / norm(gamma_new);
          
        } else {
          
          // Rotate gamma_old:
          gamma_new = abs(cos(theta) * gamma_old + sin(theta) * v * tau_root);
          gamma_new = gamma_new * tau_root / norm(gamma_new);
          
        }
        
        // Check if gamma_new is close to gamma_old:
        if(abs(gamma_new - gamma_old).max() < eps){
          
          continueLoop = false;
          
        }
        
        // Increment counter:
        j++;
        
      }
      
      // Set the solution for the current tau:
      GAMMA.col(t) = gamma_new;
      
    }
    
  }
  
  // Return output list:
  return Rcpp::List::create(Rcpp::Named("LAMBDA") = GAMMA % GAMMA);
  
}
