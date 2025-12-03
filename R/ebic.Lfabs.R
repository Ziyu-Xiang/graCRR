#---------- Initial estimator with data-driven choice of lambda2 ----------#
ebic.Lfabs <- function(X, y, G, h = 1, nlambda2 = 5, weight = NULL)
{
  x = as.matrix(X)
  if (is.null(weight))         weight = rep(1, ncol(x))
  
  # Generate the grid of candidate lambda2 values
  fit <- Lfabs(x, y, 0, G, h, weight = weight)
  lamb.max = fit$lambda[1]
  lambda2.max = 5*lamb.max
  lambda2.min = lamb.max
  lambda2.grid = exp(seq(log(lambda2.min), log(lambda2.max), length.out = nlambda2))
  
  ebic_list <- rep(0,nlambda2)
  for (llambda2 in 1:nlambda2) {
    lambda2 = lambda2.grid[llambda2]
    fit_l <- Lfabs(x, y, lambda2, G, h, weight = weight)
    opt_l <- fit_l$opt
    ebic_list[llambda2] <- fit_l$ebic[opt_l]
  }
  opt <- which.min(ebic_list)
  
  # Fit the model using the optimal lambda2
  fit = Lfabs(x, y, lambda2.grid[opt], G, h, weight = weight)
  
  val = list(Beta        = fit$Beta,
             beta        = fit$beta,
             lambda      = fit$lambda,
             lambda2     = lambda2.grid[opt],
             direction   = fit$direction,
             iter        = fit$iter,
             ebic        = fit$ebic,
             opt         = fit$opt)
  
  return(val)
}