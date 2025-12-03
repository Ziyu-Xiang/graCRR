#include <R.h>
#include <Rinternals.h>
#include <math.h>
#include <Rmath.h>
#include <stdio.h>
#include <stdlib.h>

// Define the function L_prime
double L_prime(double u) {
  if (u >= -1.0 && u <= 1.0) {
    return 1.5 * u - 0.5 * pow(u, 3.0);
  } else if (u > 1.0) {
    return 1.0;
  } else {
    return -1.0;
  }
}

void _bootstrap_individual(double *boot, double *ans, double *tmp,
  int B, double *epsilon_hat, double *X, double *W_hat, 
  double *S_inv_sqrt, double *Z, int *param){

    int n = param[0];
    int p = param[1];
    int i, j, k, b;
    double temp, sum;
    double norm_factor = (double) n * (double) (n - 1);

    for (b = 0; b < B; ++b) {

      for (k = 0; k < p; ++k) ans[k] = 0.0;

      for (i = 0; i < n; ++i) {
        for (j = 0; j < n; ++j) {
          if (i != j) {
            temp = L_prime(epsilon_hat[i] - epsilon_hat[j]) * (Z[b*n + i] + Z[b*n + j]);
            for (k = 0; k < p; ++k) {
              ans[k] += (X[i*p + k] - X[j*p + k]) * temp;
            }
          }
        }
      }

      // tmp = W_hat * ans
      for (i = 0; i < p; ++i) {
        sum = 0.0;
        for (j = 0; j < p; ++j) {
          sum += W_hat[i*p + j] * ans[j];
        }
        tmp[i] = sum;
      }

      // ans = S_inv_sqrt * tmp / norm_factor
      for (i = 0; i < p; ++i) {
        ans[i] = S_inv_sqrt[i] * tmp[i] / norm_factor;
      }

      // Store |ans_j| into boot (p x B matrix)
      for (j = 0; j < p; ++j) {
        boot[j + p * b] = fabs(ans[j]);
      }
    }
}

SEXP bootstrap_individual(SEXP B_, SEXP epsilon_hat_, SEXP X_, SEXP W_hat_,
                          SEXP S_hat_sqrt_inverse_, SEXP Z_, SEXP Param) {
  int B, p;
  int *param;
  double *epsilon_hat, *X, *W_hat, *S_inv_sqrt, *Z;
  
  B           = INTEGER(B_)[0];
  param       = INTEGER(Param);
  p           = param[1];
  epsilon_hat = REAL(epsilon_hat_);
  X           = REAL(X_);
  W_hat       = REAL(W_hat_);
  S_inv_sqrt  = REAL(S_hat_sqrt_inverse_);
  Z           = REAL(Z_);

  /* allocate result matrix p x B */
  SEXP bootstraps = PROTECT(allocMatrix(REALSXP, p, B));
  double *boot = REAL(bootstraps);

  double *ans, *tmp;
  ans  = (double*)calloc(p, sizeof(double));
  tmp  = (double*)calloc(p, sizeof(double));

  _bootstrap_individual(boot, ans, tmp, B, epsilon_hat, X, W_hat, S_inv_sqrt, Z, param);

  free(ans);
  free(tmp);

  UNPROTECT(1);
  return bootstraps;
}

