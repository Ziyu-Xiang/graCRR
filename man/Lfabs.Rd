\name{Lfabs}
\alias{Lfabs}
\title{
    A forward and backward stagewise algorithm for high-dimensional convoluted rank regression with sparse Laplacian shrinkage
}
\description{
    Fit a high-dimensional linear model using the convoluted rank loss with
    \eqn{\ell_1} and graph Laplacian penalties. The function computes a solution path
    over the \eqn{\ell_1} tuning parameter \code{lambda} while fixing the Laplacian penalty
    parameter \code{lambda2}, and selects the optimal model along the path using EBIC.
    The graph structure is encoded by \code{G} and determines which predictors are
    penalized jointly.
}
\usage{
Lfabs(X, y, lambda2, G, h = 1, stopping = TRUE, eps = 0.02, xi = 1e-6,
    iter = 1e4, lambda.min = 1e-4, weight = NULL)
}
\arguments{
  \item{X}{
    An \eqn{n \times p} design matrix of predictors.
  }
  \item{y}{
    A numeric response vector of length \eqn{n}.
  }
  \item{lambda2}{
    Nonnegative tuning parameter for the graph Laplacian penalty term
    \eqn{\frac{1}{2}\lambda_2 \boldsymbol{\beta}^\top L \boldsymbol{\beta}},
    where \eqn{L} is the graph Laplacian matrix derived from \code{G}. 
    This value is held fixed while the \eqn{\ell_1} tuning parameter varies along a path.
  }
  \item{G}{
	A square matrix representing the weighted graph structure of the predictors.
    Each entry \eqn{G[i, j]} indicates the weight of the edge between node \eqn{i} and node \eqn{j}.
    A value of zero means no edge is present. Diagonal entries are ignored.
  }
  \item{h}{
    A positive bandwidth parameter used in the smoothed (convoluted) rank loss
    and its derivatives. Default is \code{1}.
  }
  \item{stopping}{
    Logical indicator of whether to stop the iterations when \code{lambda} is less than
    \code{lambda.min}. Default is \code{TRUE}.
  }
  \item{eps}{
    The step size for updating coefficients. Default is \code{0.02}.
  }
  \item{xi}{
    A small positive threshold used in the algorithm. Default is \code{1e-6}.
  }
  \item{iter}{
    The maximum number of iterations allowed. Default is \code{1e4}.
  }
  \item{lambda.min}{
    The smallest value for \code{lambda}, used as a stopping criterion for the solution path.
    Default is \code{1e-4}.
  }
  \item{weight}{
    A weight vector of length \eqn{p} for the adaptive lasso penalty. Default is \code{NULL}, 
    which corresponds to weight 1 for all \eqn{\beta_j}, i.e., the standard lasso.
  }
}
\value{
A list with the following components:
\item{Beta}{
    A matrix containing the coefficient estimates along the solution path, one column for each step.
}
\item{beta}{
    A numeric vector of length \eqn{p} giving the coefficient estimates at the
    EBIC-selected optimal tuning parameter.
}
\item{lambda}{
    A numeric vector containing the sequence of \eqn{\ell_1} tuning parameters
    along the path.
}
\item{direction}{
    A numeric vector indicating the sign of the change in the updated
    coefficient at each iteration. \code{1} indicates a forward step; 
    \code{0} indicates a backward step.
}
\item{iter}{
    Integer scalar giving the number of iterations.
}
\item{ebic}{
    A numeric vector of EBIC values along the solution path.
}
\item{opt}{
    Integer index of the EBIC-selected optimal model along the path.
}
}
\details{
The Laplacian shrinkage level is fixed, which means \code{lambda2} is predetermined.
EBIC is evaluated along the path and is used to select the optimal model.
}
\references{
Shi, X., Huang, Y., Huang, J., and Ma, S. (2018).
A forward and backward stagewise algorithm for nonconvex loss functions with adaptive lasso.
\emph{Computational Statistics & Data Analysis}, \bold{124}, 235--251.

Chen, J. and Chen, Z. (2008).
Extended Bayesian information criteria for model selection with large model spaces.
\emph{Biometrika}, \bold{95}(3), 759--771.

Tan, X., Zhang, X., Cui, Y., and Liu, X. (2024).
Uncertainty quantification in high-dimensional linear models incorporating graphical structures with applications to gene set analysis.
\emph{Bioinformatics}, \bold{40}(9), btae541.
}
\examples{
set.seed(123)

## Simulate data (see Data_simu for details)
dat <- Data_simu(n = 50, p = 100, s0 = 5,
                 xcov = "AR1", error = "norm", signal = 1)
X <- dat$X
y <- dat$y
G <- dat$G

fit <- Lfabs(X, y, lambda2 = 1, G = G, h = 1)
beta_hat <- fit$beta
}