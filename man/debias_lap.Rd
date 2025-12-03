\name{debias_lap}
\alias{debias_lap}
\title{Graph-assisted debiased estimator for high-dimensional linear models}
\description{
Compute a debiased estimator for high-dimensional linear models fitted
using the convoluted rank loss with \eqn{\ell_1} and Laplacian penalties.
}
\usage{
	debias_lap(X, y, beta, G, h)
}
\arguments{
	\item{X}{
    An \eqn{n \times p} design matrix of predictors.
  }
    \item{y}{
    A numeric response vector of length \eqn{n}.
  }
    \item{beta}{
    A numeric vector of length \eqn{p} giving an initial estimator of the regression
    coefficients, i.e., the solution to a graph-assisted convoluted rank
    regression problem with \eqn{\ell_1} and Laplacian penalties.
  }
	\item{G}{
	A square matrix representing the weighted graph structure of the predictors.
    Each entry \eqn{G[i, j]} indicates the weight of the edge between node \eqn{i} and node \eqn{j}.
    A value of zero means no edge is present. Diagonal entries are ignored.
  }
	\item{h}{
    A positive bandwidth parameter used in the smoothed (convoluted) rank loss
    and its derivatives. It should match the bandwidth used to obtain
    the initial estimator \code{beta}.
  }
}
\value{
	A numeric vector of length \eqn{p} giving the debiased estimator of the regression
coefficients.
}
\references{
Tan, X., Zhang, X., Cui, Y., and Liu, X. (2024).
Uncertainty quantification in high-dimensional linear models incorporating graphical structures with applications to gene set analysis.
\emph{Bioinformatics}, \bold{40}(9), btae541.

Zhang, C.-H. and Zhang, S. S. (2014).
Confidence intervals for low dimensional parameters in high dimensional linear models.
\emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)},
\bold{76}(1), 217--242.
}
\examples{
set.seed(123)

## Simulate data (see Data_simu for details)
dat <- Data_simu(n = 50, p = 100, s0 = 5,
                 xcov = "AR1", error = "norm", signal = 1)
X <- dat$X
y <- dat$y
G <- dat$G

## Suppose we start from a simple initial estimator
beta_init <- rep(0, ncol(X))

beta_db <- debias_lap(X, y, beta_init, G, h = 1)
}