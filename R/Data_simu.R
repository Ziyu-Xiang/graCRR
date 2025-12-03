#---------- Simulate data ----------#
Data_simu <- function(n, p, s0, xcov, error, signal=1) {
  # Check if p is a multiple of 5
  if (p %% 5 != 0) {
    stop(paste("The program stops running because", p, "is not a multiple of 5."))
  }

  Covariance_AR <- function(p, rho) {
    cluster_size <- 5
    K <- p / cluster_size
    
    # AR(1) covariance for a single cluster
    C <- matrix(0, cluster_size, cluster_size)
    for (i in 1:(cluster_size-1)) {
      for (j in (i + 1):cluster_size) {
        C[i, j] <- rho^(abs(i - j))
      }
    }
    
    Cov <- C + t(C)
    diag(Cov) <- rep(1, cluster_size)
    
    # Block-diagonal covariance with K clusters
    Sigma <- Matrix::bdiag(replicate(K, Cov, simplify = FALSE))
    Sigma <- as.matrix(Sigma)
    
    return(Sigma)
  }
  
  Covariance_band <- function(p, rho_seq) {
    # Normalize rho_seq to have unit L2 norm
    rho_seq = rho_seq/sqrt(sum(rho_seq^2))
    n.rho = length(rho_seq)
    K <- p / n.rho
    
    Cov <- matrix(0, n.rho, n.rho)
    for (i in 1:(n.rho-1)) {
      tmp = rho_seq[(1+i):n.rho]%*% rho_seq[1:(n.rho-i)]
      Cov[row(Cov) == col(Cov) + i] = tmp
      Cov[row(Cov) + i == col(Cov)] = tmp
    }
    diag(Cov) <- rep(1, n.rho)
    
    # Block-diagonal covariance with K clusters
    Sigma <- Matrix::bdiag(replicate(K, Cov, simplify = FALSE))
    Sigma <- as.matrix(Sigma)
    
    return(Sigma)
  }
  
  #---------- Generate Data ----------#
  Data <- function(n, p, s0, Cov, epsilon, signal) {
    beta <- rep(0, p)
    beta[1:s0] <- signal
    X <- rmnorm(n, rep(0, p), Cov)
    y = X %*% beta + epsilon
    list(X = X, y = y, beta = beta)
  }
  
  # Generate error term
  if (error == "norm") {
    # Standard normal errors
    epsilon <- rnorm(n, mean = 0, sd = 1)
  } else if (error == "mnorm") {
    # Mixture normal errors: 95% from N(0, 1), 5% from N(0, 100)
    z <- rbinom(n, 1, prob = 0.95)
    epsilon <- rnorm(n, 0, 1) * z + rnorm(n, 0, 100) * (1 - z)
  } else if (error == "t3") {
    # t-distributed errors with 3 degrees of freedom (heavy tails)
    epsilon <- rt(n, df = 3)
  } else if (error == "cauchy") {
    # t-distributed errors with 1 degree of freedom (= standard Cauchy)
    epsilon <- rt(n, df = 1)
  } else {
    stop("No valid error specification.")
  }
  
  # Generate design covariance matrix
  if (xcov == "AR1") {
    # AR(1) covariance with rho = 0.5
    rho = 0.5
    Cov = Covariance_AR(p, rho)
  } else if (xcov == "AR2") {
    # AR(1) covariance with rho = 1 / 1.01
    rho = 1/1.01
    Cov = Covariance_AR(p, rho)
  } else if (xcov == "Band") {
    # Banded (moving average–type) covariance structure
    rho_seq = runif(5)
    Cov = Covariance_band(p, rho_seq)
  } else {
    stop("No valid covariance specification.")
  }
  
  Theta = solve(Cov)
  G = ifelse(abs(Theta)>1e-8,1,0) # Theoretically symmetric
  
  Da = Data(n, p, s0, Cov, epsilon, signal)
  
  list(X = Da$X, y = Da$y, beta_T=Da$beta, G = G)
}
