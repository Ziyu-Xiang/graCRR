#---------- Inference ----------#
inference = function(X, y, beta_hat, beta_true, debias, G, h, alpha=0.05, B=1000){
  
  pairs_full <- function(n) {
    # fast construction: create all ordered pairs (i,j) with i != j
    # returns 2 x (n*(n-1)) matrix: first row i, second row j
    ii <- rep(1:n, each = n)
    jj <- rep(1:n, times = n)
    keep <- ii != jj
    mat <- cbind(ii[keep], jj[keep])
    # transpose to 2 x m
    t(mat)
  }
  
  n = nrow(X)
  p = ncol(X)
  epsilon_hat = y-X%*%beta_hat
  
  pairs <- pairs_full(n)
  delta_y <- y[pairs[1, ]] - y[pairs[2, ]]
  delta_x <- X[pairs[1, ], ] - X[pairs[2, ], ]
  # Residual differences u for pairwise observations
  u <- delta_y - as.vector(delta_x %*% beta_hat)
  
  inv_gram <- calculate_inv_gram(X,G)$inverse_Gram
  term2 <- 2*sum(L_h_double_prime(u, h))
  W_hat <-  1/term2 * t(inv_gram) * (n*(n-1))
  
  S_hat_sqrt_diag  =  sqrt(diag(W_hat))
  S_hat_sqrt_inverse = matrix(0,p,p)
  diag(S_hat_sqrt_inverse) = 1/sqrt(diag(W_hat))  
  
  bootstraps = bootstrap_individual(B=B,epsilon_hat,X,W_hat,S_hat_sqrt_inverse)
  Q = apply(bootstraps, 1, quantile, 1 - alpha)
  
  val = list(cover = (debias - S_hat_sqrt_diag * Q < beta_true) & 
               (debias + S_hat_sqrt_diag * Q > beta_true),
             length = 2*S_hat_sqrt_diag * Q)
  
  return(val)
}
