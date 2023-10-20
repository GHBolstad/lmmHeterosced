#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{

  // Data
  DATA_MATRIX(Y);        // Response variable
  DATA_SPARSE_MATRIX(X); // Fixed effects design matrix
  DATA_SPARSE_MATRIX(Z); // Random effects design matrix
  DATA_INTEGER(n_random);
  DATA_FACTOR(u_start);
  DATA_FACTOR(u_end);
  DATA_MATRIX(Q);        // Design matrix for the log residual variance
  DATA_VECTOR(SE);       // Standard error of each element of Y


  // Parameters
  PARAMETER_MATRIX(b);        // The fixed effects
  PARAMETER_MATRIX(b_ln_R);   // The linear parameters the log residual variance
  PARAMETER_MATRIX(u);        // The random effects
  PARAMETER_VECTOR(ln_G);     // The log variances of the random effects

  Type nll=0;
  int n = Y.rows();

  // Mixed model equation
  vector<Type> residuals = Y - X*b - Z*u;

  //----------
  // Residuals
  //----------
  // log residual variance
  vector<Type> ln_R = Q*b_ln_R;
  
  int i;
  for(i=0; i<n; i++){
    nll -= dnorm(residuals(i), Type(0.0), exp(Type(0.5)*ln_R(i))+SE(i), true);
  }


  //---------------
  // Random effects
  //---------------

  int k;
  for(k=0; k<n_random; k++){
    for(i=u_start(k); i<u_end(k); i++){
      nll -= dnorm(u(i), Type(0.0), exp(Type(0.5)*ln_G(k)), true);
    }
  }
  
  //-------
  // Report
  //-------

  REPORT(residuals);

  return nll;
}



