#include <Rcpp.h>
#include <RcppEigen.h>
using namespace Rcpp;
using namespace Eigen;

// Define the function L_prime
double L_prime(double u) {
  if (u >= -1.0 && u <= 1.0) {
    return 1.5 * u - 0.5 * std::pow(u, 3.0);
  } else if (u > 1.0) {
    return 1.0;
  } else {
    return -1.0;
  }
}

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
NumericMatrix bootstrap_individual(
  const int B,
  const Eigen::VectorXd& epsilon_hat,
  const Eigen::MatrixXd& X,
  const Eigen::MatrixXd& W_hat,
  const Eigen::MatrixXd& S_hat_sqrt_inverse
) {
  const int n = X.rows();
  const int p = X.cols();
  
  // Each column is one bootstrap draw; each row corresponds to a coefficient
  NumericMatrix bootstraps(p, B);
  
  for (int b = 0; b < B; ++b) {
    // Generate multiplier weights Z
    NumericVector z = Rcpp::rnorm(n);
    Map<VectorXd> Z(z.begin(), n);
    Eigen::VectorXd ans = Eigen::VectorXd::Zero(p);
    
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        if (i != j) {
          double lij = L_prime(epsilon_hat(i) - epsilon_hat(j));
          double zij = Z(i) + Z(j);
          double s   = lij * zij; // scalar

          // (X.row(i) - X.row(j)) is a RowVector
          RowVectorXd diff = X.row(i) - X.row(j); // 1×p
          ans.noalias() += diff.transpose() * s; // p×1
        }
      }
    }

    ans = (S_hat_sqrt_inverse * W_hat * ans) / double(n * (n - 1));

    for (int j = 0; j < p; ++j) {
      // bootstraps(j, b) = ans(j);
      bootstraps(j, b) = std::abs(ans(j));
    }
  }

  return bootstraps;
}

