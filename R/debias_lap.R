#---------- Debias ----------#
debias_lap <- function(X, y, beta, G, h){
  
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
  pairs <- pairs_full(n)
  delta_y <- y[pairs[1, ]] - y[pairs[2, ]]
  delta_x <- X[pairs[1, ], ] - X[pairs[2, ], ]
  # Residual differences u for pairwise observations
  u <- delta_y - as.vector(delta_x %*% beta)
  
  inv_gram <- calculate_inv_gram(X,G)$inverse_Gram
  
  weights <- L_h_prime(u, h)
  weighted_delta_x <- sweep(delta_x, 1, weights, FUN = "*")
  term1 <- colSums(weighted_delta_x)
  term2 <- 2*sum(L_h_double_prime(u, h))
  debias_Lap <- beta + 1/term2 * t(inv_gram) %*% term1
  
  return(debias_Lap)
}
