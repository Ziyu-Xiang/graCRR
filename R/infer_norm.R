#---------- Inference ----------#
infer_norm = function(X, y, beta_hat, beta_true = NULL, debias, G, h, 
                      alpha=0.05, inv_method = "Tan2024"){
  
  n <- nrow(X)
  Lf <- make_Lh_functions(h)
  if (inv_method == "Tan2024"){
    Theta <- calculate_inv_gram(X, G)$inverse_Gram
  }
  else if (inv_method == "Cai2025"){
    Theta <- flare::sugm(X, verbose = F, nlambda = 2, lambda.min.ratio = 0.25)
    Theta <- Theta$icov[[Theta$nlambda]]
  }
  else {
    stop("No valid method specification.")
  }
  theta_jj <- diag(Theta)
  
  eps <- as.numeric(y - X %*% beta_hat)
  factor <- multi_factor(eps, Lf$Lh1, Lf$Lh2)
  var_Lap <- factor * theta_jj
  sehat <- sqrt(var_Lap/n)
  Delta <- qnorm(1-alpha/2)*sehat
  
  p_value <- 2 * pnorm(abs(debias / sehat), lower.tail = FALSE)
  
  if (is.null(beta_true)) {
    val = list(
      cover    = NULL,
      length   = 2*Delta,
      lower    = debias-Delta,
      upper    = debias+Delta,
      Inv_Gram = Theta,
      sehat    = sehat,
      p_value  = p_value
    )
  } else {
    val = list(
      cover    = (debias-Delta<=beta_true) & 
        (debias+Delta>=beta_true),
      length   = 2*Delta,
      lower    = debias-Delta,
      upper    = debias+Delta,
      Inv_Gram = Theta,
      sehat    = sehat,
      p_value  = p_value
    )
  }
  
  return(val)
}
