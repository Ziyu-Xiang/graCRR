\name{ebic.Lfabs}
\alias{ebic.Lfabs}
\title{
	A forward and backward stagewise algorithm for high-dimensional convoluted rank regression with sparse Laplacian shrinkage
}
\description{
    Fit a high-dimensional linear model using the convoluted rank loss with
    \eqn{\ell_1} and graph Laplacian penalties, where the Laplacian penalty
    parameter \code{lambda2} is selected in a data-driven way using EBIC.
    The function evaluates a grid of candidate \code{lambda2} values, applies
    \code{Lfabs} at each value, and chooses the Laplacian penalty that yields
    the smallest EBIC. The graph structure is encoded by \code{G}.
}
\usage{
	ebic.Lfabs(X, y, G, h = 1, nlambda2 = 5, weight = NULL)
}
\arguments{
  \item{X}{
    An \eqn{n \times p} design matrix of predictors.
  }
  \item{y}{
    A numeric response vector of length \eqn{n}.
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
  \item{nlambda2}{
	The number of candidate \code{lambda2} values in the grid for EBIC-based
    selection. Default is \code{5}.
  }
  \item{weight}{
    A weight vector of length \eqn{p} for the adaptive lasso penalty. Default is \code{NULL}, 
    which corresponds to weight 1 for all \eqn{\beta_j}, i.e., the standard lasso.
  }
}
\value{
A list with the following components:
\item{Beta}{
    A matrix containing the coefficient estimates along the solution path
	(for the selected \code{lambda2}), one column for each step.
}
\item{beta}{
    A numeric vector of length \eqn{p} giving the coefficient estimates at the
    EBIC-selected optimal \eqn{\ell_1} tuning parameter for the chosen \code{lambda2}.
}
\item{lambda}{
    A numeric vector containing the sequence of \eqn{\ell_1} tuning parameters
    along the path for the selected \code{lambda2}.
}
\item{lambda2}{
    The selected Laplacian penalty parameter \code{lambda2} that minimizes EBIC
    over the candidate grid.
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
    A numeric vector of EBIC values along the solution path (for the selected
    \code{lambda2}).
}
\item{opt}{
    Integer index of the EBIC-selected optimal model along the path.
}	
}
\details{
The function first fits the model with \code{lambda2 = 0} via \code{Lfabs} to
obtain a reference \eqn{\ell_1} tuning parameter. It then constructs a grid of
candidate Laplacian penalties \code{lambda2} on a logarithmic scale and, for each
candidate, calls \code{Lfabs} and records the minimum EBIC along the corresponding
solution path. The \code{lambda2} value with the smallest EBIC is selected, and
the final model is refitted at this optimal Laplacian penalty.
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

fit <- ebic.Lfabs(X, y, G, h = 1, nlambda2 = 5)
beta_hat <- fit$beta
lambda2_hat <- fit$lambda2
}
\seealso{\code{\link{Lfabs}}}