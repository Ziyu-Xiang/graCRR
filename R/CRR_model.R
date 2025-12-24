L_h <- function(u, h) {
  result <- numeric(length(u))
  result[u >= h] <- u[u >= h]
  result[u <= -h] <- -u[u <= -h]
  idx <- which(abs(u) < h)
  result[idx] <- (3 * u[idx]^2) / (4 * h) - (u[idx]^4) / (8 * h^3) + (3 * h) / 8
  
  return(result)
}

L_h_prime <- function(u, h) {
  grad <- numeric(length(u))
  grad[u >= h] <- 1
  grad[u <= -h] <- -1
  idx <- which(abs(u) < h)
  grad[idx] <- (3 * u[idx]) / (2 * h) - (u[idx]^3) / (2 * h^3)
  
  return(grad)
}

L_h_double_prime <- function(u, h) {
  second_derivative <- numeric(length(u)) 
  
  idx <- which(abs(u) < h)
  second_derivative[idx] <- 3 / (2 * h) - (3 * u[idx]^2) / (2 * h^3)
  
  return(second_derivative)
}

make_Lh_functions <- function(h) {
  
  Lh1 <- function(u) {
    matrix(
      L_h_prime(as.vector(u), h),
      nrow = nrow(u),
      ncol = ncol(u)
    )
  }
  
  Lh2 <- function(u) {
    matrix(
      L_h_double_prime(as.vector(u), h),
      nrow = nrow(u),
      ncol = ncol(u)
    )
  }
  
  list(Lh1 = Lh1, Lh2 = Lh2)
}

multi_factor <- function(eps, Lh1, Lh2) {
  eps <- as.numeric(eps)
  n <- length(eps)
  if (n < 2) stop("Need at least n>=2 residuals.")
  
  D <- outer(eps, eps, "-")     # D[i,j] = eps[i] - eps[j]
  off <- !diag(n)              # i != j
  
  # B_hat = average of L_h'' over i != j
  L2 <- Lh2(D)
  B_hat <- mean(L2[off])
  
  # m_hat[i] = average of L_h' over j != i, conditional on eps_i
  L1 <- Lh1(D)
  m_hat <- rowSums(L1 * off) / (n - 1)
  
  # A_hat = average of m_hat[i]^2
  A_hat <- mean(m_hat^2)
  
  if (abs(B_hat) < .Machine$double.eps) stop("B_hat is numerically zero; variance blows up.")
  factor <- (A_hat / (B_hat^2))
  
  return(factor)
}

calculate_inv_gram <- function(x, G, k = NULL) {
  n <- nrow(x)
  p <- ncol(x)
  if (is.null(k)) k <- p
  Z <- matrix(0, n, p)
  diag(G) <- 0
  inverse_Gram <- matrix(0, p, p)
  for (j in 1:p) {
    D <- which(G[j, ] != 0)
    if (length(D) == 0) {
      inverse_Gram[, j] <- 0
      Z[, j] <- x[, j]
      inverse_Gram[j, j] <- n / sum(x[, j]^2)
    } else {
      if (length(D) >= 1 & length(D) <= k) {
        fit <- lm(x[, j] ~ x[, D] - 1)
        gamma <- fit$coefficients
        Z[, j] <- fit$residuals
        tao <- sum((fit$residuals)^2) / n
        inverse_Gram[j, j] <- 1 / tao
        inverse_Gram[, j][D] <- -gamma / tao
      }
      
      if (length(D) > k) {
        glmnetfit <- cv.glmnet(x[, D], x[, j])
        gamma <- as.vector(predict(glmnetfit, x[, D],
                                   type = "coefficients",
                                   s = "lambda.1se"
        ))[-1]
        prediction <- predict(glmnetfit, x[, D], s = "lambda.1se")
        Z[, j] <- x[, j] - prediction
        tao <- as.numeric((x[, j] %*%
                             (x[, j] - predict(glmnetfit, x[, D], s = "lambda.1se"))) / n)
        inverse_Gram[j, j] <- 1 / tao
        inverse_Gram[, j][D] <- -gamma / tao
      }
    }
  }
  list(inverse_Gram = inverse_Gram, Z = Z)
}
