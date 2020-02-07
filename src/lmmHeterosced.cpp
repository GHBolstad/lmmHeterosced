#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{

  // Data
  DATA_MATRIX(Y);        // The response varaible
  DATA_SPARSE_MATRIX(X); // The fixed effects design matrix
  DATA_SPARSE_MATRIX(Z); // The random effects design matrix
  DATA_INTEGER(n_random);
  DATA_FACTOR(u_start);
  DATA_FACTOR(u_end);
  DATA_MATRIX(Q);        // Design matrix for the residual variance


  // Parameters
  PARAMETER_MATRIX(b);        // The fixed effects
  PARAMETER_MATRIX(u);        // The random effects
  PARAMETER_MATRIX(b_ln_R);   // The linear parameters the log residual variance
  PARAMETER_VECTOR(ln_G);     // log variances of the random effects


  int i, k;
  Type nll=0;
  int n = Y.rows();

  // Mixed model equation
  vector<Type> residuals = Y - X*b - Z*u;


  //----------
  // Residuals
  //----------
  // log residual variance
  vector<Type> ln_R = Q*b_ln_R;

  for(i=0;i<n;i++){
    nll -= dnorm(residuals(i), Type(0.0), exp(Type(0.5)*ln_R(i)), true);
  }


  //---------------
  // Random effects
  //---------------

  for(k=0; k<n_random; k++){
    for(i=u_start(k);i<u_end(k);i++){
      nll -= dnorm(u(i), Type(0.0), exp(Type(0.5)*ln_G(k)), true);
    }
  }

  ADREPORT(residuals);

  return nll;
}



