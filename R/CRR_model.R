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