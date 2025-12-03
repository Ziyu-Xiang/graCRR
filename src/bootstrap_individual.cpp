#include <Rcpp.h>
#include <RcppEigen.h>
using namespace Rcpp;
using namespace Eigen;

// Define the function L_prime
double L_prime(double u) {
  if(u >= -1 && u <= 1) {
    return 3.0/2.0*u - pow(u, 3)/2;
  } else if (u > 1) {
    return 1.0;
  } else {
    return -1.0;
  }
}

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
NumericMatrix bootstrap_individual(int B, Eigen::VectorXd epsilon_hat, Eigen::MatrixXd X, Eigen::MatrixXd W_hat, Eigen::MatrixXd S_hat_sqrt_inverse) {
  int n = X.rows();
  int p = X.cols();
  
  NumericMatrix bootstraps(p, B); // ← 每列是一次bootstrap，每行为某个系数
  
  for (int b = 0; b < B; ++b) {
    NumericVector z = rnorm(n);
    Map<VectorXd> Z(as<Map<VectorXd>>(z));
    Eigen::VectorXd ans = Eigen::VectorXd::Zero(p);
    
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        if (i != j) {
          ans += L_prime(epsilon_hat[i] - epsilon_hat[j]) * (X.row(i) - X.row(j)) * (Z[i] + Z[j]);
        }
      }
    }

    ans = S_hat_sqrt_inverse * W_hat * ans / (n * (n - 1));

    for (int j = 0; j < p; ++j) {
      // bootstraps(j, b) = ans(j);  // 存储第 j 个系数的偏差项
      bootstraps(j, b) = std::abs(ans(j));
    }
  }

  return bootstraps;
}

