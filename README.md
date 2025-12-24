# graCRR
**graCRR** is an R package for **graph-assisted convoluted rank regression** in high-dimensional linear models.

It provides tools for:

- forward–backward stagewise (Fabs) estimation with lasso and Laplacian penalties
- EBIC-based selection of the Laplacian penalty
- debiasing for inference
- inference based on asymptotic normality, including standard errors, confidence intervals, and p-values
- optional multiplier bootstrap procedures for interval construction
- simulation tools for method evaluation

The proposed methods are robust to heavy-tailed errors and incorporate structural information through graphical relationships among predictors.

# Installation

    #install Rtools 4.4 (https://cran.r-project.org/bin/windows/Rtools)
    #install.packages("devtools")
    #install.packages("Rcpp")
    library(devtools)
    install_github("Ziyu-Xiang/graCRR")

# Usage

- [x] [graCRR-manual](https://github.com/Ziyu-Xiang/graCRR/tree/main/inst/graCRR-manual.pdf) ------------ Details of the usage of the package.

# Example

    library(graCRR)
    
    set.seed(123)
    
    ## Simulate data
    dat <- Data_simu(n = 80, p = 100, s0 = 5,
                     xcov = "AR1", error = "norm", signal = 1)
    X <- dat$X
    y <- dat$y
    G <- dat$G
    beta_true <- dat$beta_T
    
    ## Step 1: EBIC selection of Laplacian penalty
    fit <- ebic.Lfabs(X, y, G, h = 1)
    beta_hat <- fit$beta
    
    ## Step 2: Debiasing
    beta_db <- debias_lap(X, y, beta_hat, G, h = 1)
    
    ## Step 3: Inference based on normal approximation
    res_norm <- infer_norm(X, y, beta_hat,
                            beta_true = beta_true,
                            debias = beta_db,
                            G = G, h = 1,
                            inv_method = "Tan2024")
    
    res_norm$cover    # coverage indicators (simulation only)
    res_norm$length   # confidence interval lengths
    res_norm$p_value  # two-sided p-values for H0: beta_j = 0


# References
- Tan, X., Zhang, X., Cui, Y., & Liu, X. (2024).
   *Uncertainty quantification in high-dimensional linear models incorporating graphical structures with applications to gene set analysis.*
   **Bioinformatics**, **40**(9), btae541.

- Cai, L., Guo, X., Lian, H., & Zhu, L. (2025).
   *Statistical inference for high-dimensional convoluted rank regression.*
   **Journal of the American Statistical Association**, 1–25.

- Shi, X., Huang, Y., Huang, J., & Ma, S. (2018).
   *A forward and backward stagewise algorithm for nonconvex loss functions with adaptive lasso.*
   **Computational Statistics & Data Analysis**, **124**, 235–251.

- Zhang, C.-H., & Zhang, S. S. (2014).
  *Confidence intervals for low dimensional parameters in high dimensional linear models.*
  **Journal of the Royal Statistical Society: Series B (Statistical Methodology)**, **76**(1), 217–242.

- Chen, J., & Chen, Z. (2008).
   *Extended Bayesian information criteria for model selection with large model spaces.*
   **Biometrika**, **95**(3), 759–771.

# Development
The R-package is developed by Ziyu Xiang (xiangziyu@stu.sufe.edu.cn), Xiangyong Tan, Shishuo Guo, Xu Liu.



