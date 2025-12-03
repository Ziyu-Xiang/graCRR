#---------- Initial estimator with a specified lambda2 ----------#
Lfabs <- function(X, y, lambda2, G, h = 1,
                  stopping = TRUE, eps = 0.02, xi = 1e-6,
                 iter = 1e4, lambda.min = 1e-4, weight = NULL) {
  
  Adj <- function(G) {
    p <- ncol(G)
    A <- matrix(0, p, p)
    A <- ifelse(G != 0, 1, 0)
    diag(A) <- 0
    
    degree = colSums(abs(A))
    adj = c(A[upper.tri(A)])
    id = which(adj!=0)
    adj = adj[id]
    adjRow = c(c(sequence(sequence(p-1)))-1)[id]
    adjCol = c(rep(1:(p-1), 1:(p-1)))[id]
    
    val = list(edge    = adj,
               edgerow = adjRow,
               edgecol = adjCol,
               degree  = degree,
               ledge   = length(adj))
    
    return(val)
  }
  
  x = as.matrix(X)
  n = nrow(x)
  p = ncol(x)
  max_s = as.integer(n/log(n))
  
  # Standardize X and center y
  x <- scale(x, center = TRUE, scale = TRUE)
  sdx   <- attr(x, "scaled:scale")
  y <- scale(y, center = TRUE, scale = FALSE)
  
  if (is.null(lambda.min)) lambda.min = {if (n > p) 1e-4 else 0.02}
  if (is.null(weight))     weight = rep(1, p)
  
  # Adjacency information derived from G
  adj = Adj(G)
  param = c(n, p, adj$ledge, 0, 0, 0, 1)
  
  ebic.fit <- .Call("LFabs",
                    as.numeric(y),
                    as.numeric(t(x)),
                    as.numeric(weight),
                    as.numeric(eps),
                    as.numeric(lambda.min),
                    as.numeric(xi),
                    as.integer(stopping),
                    as.integer(iter),
                    as.integer(param),
                    as.numeric(adj$edge),
                    as.integer(adj$edgerow),
                    as.integer(adj$edgecol),
                    as.numeric(adj$degree),
                    as.integer(max_s),
                    as.numeric(lambda2),
                    as.numeric(h) )
  
  fit <- ebic.fit[[1]]
  fit$opt <- which.min(fit$bic) # index of the optimal EBIC along the path
  theta <- Matrix::sparseMatrix(fit$index_i, fit$index_j,
                                x = fit$beta, dims = c(p, fit$iter), index1 = FALSE)
  
  val <- list(
    Beta      = theta/sdx, # matrix of beta estimates at each iteration
    beta      = theta[, fit$opt]/sdx, # beta at the optimal lambda (selected by EBIC)
    lambda    = fit$lambda, # lambda path (vector)
    direction = fit$direction, # sign of change in beta_k at each iteration (vector)
    iter      = fit$iter, # number of iterations (integer scalar)
    ebic      = fit$bic, # EBIC value at each iteration (vector)
    opt       = fit$opt # index of the optimal EBIC (integer scalar)
  )
  
  return(val)
}
