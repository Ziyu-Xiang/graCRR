\name{Data_simu}
\alias{Data_simu}
\title{Simulate high-dimensional linear model data}
\description{
Generate data from a linear regression model with a structured covariance matrix for the design,
and various choices of error distributions.
}
\usage{
	Data_simu(n, p, s0, xcov, error, signal=1)
}
\arguments{
	\item{n}{
    Sample size.
	}
	\item{p}{
    Number of predictors \eqn{p}. It must be a multiple of 5; otherwise the function stops with an error.
    }
	\item{s0}{
    Number of nonzero coefficients in the true regression vector \eqn{\boldsymbol{\beta}}.
    The first \code{s0} coordinates of \eqn{\boldsymbol{\beta}} are set to \code{signal}, and the remaining are zeros.
	}
	\item{xcov}{
    A character string specifying the covariance structure of the design matrix \eqn{\mathbf{X}}.
    Possible values are:
    \describe{
      \item{\code{"AR1"}}{Block-diagonal AR(1) covariance with cluster size 5 and correlation parameter \eqn{\rho = 0.5}.}
      \item{\code{"AR2"}}{Block-diagonal AR(1) covariance with cluster size 5 and correlation parameter \eqn{\rho = 1/1.01}, 
	  corresponding to a highly correlated design.}
      \item{\code{"Band"}}{Block-diagonal banded (moving-average-type) covariance. 
	  A random vector \eqn{\boldsymbol{\rho}} of length 5 is generated from \code{runif(5)} and normalized to have unit \eqn{\ell_2}-norm; 
	  the within-block covariance is then constructed from \eqn{\boldsymbol{\rho}} and replicated along the diagonal.}
	} 
	}
	\item{error}{
	A character string specifying the distribution of the error term \eqn{\boldsymbol{\varepsilon}}.
    Possible values are:
    \describe{
      \item{\code{"norm"}}{Standard normal errors, \eqn{\varepsilon_i \sim N(0, 1)}.}
      \item{\code{"mnorm"}}{Mixture normal errors: with probability 0.95, \eqn{N(0,1)}; with probability 0.05, \eqn{N(0,100)}.}
      \item{\code{"t3"}}{t-distributed errors with 3 degrees of freedom, \eqn{\varepsilon_i \sim t_3}.}
      \item{\code{"cauchy"}}{t-distributed errors with 1 degree of freedom, i.e., standard Cauchy errors.}
    }
	}
	\item{signal}{
    Signal strength for the nonzero coefficients. The first \code{s0} entries of the true
    regression vector \eqn{\boldsymbol{\beta}} are set to \code{signal}. Default is \code{1}.
  }
}
\value{
	A list with the following components:
  \item{X}{
    An \eqn{n \times p} design matrix simulated from a multivariate normal distribution
    with mean zero and the specified covariance matrix.
  }
  \item{y}{
    The response vector of length \eqn{n}, generated as \eqn{\mathbf{y} = \mathbf{X}\boldsymbol{\beta} + \boldsymbol{\varepsilon}}.
  }
  \item{beta_T}{
    The true regression coefficient vector \eqn{\boldsymbol{\beta}} of length \eqn{p},
    where the first \code{s0} entries equal \code{signal} and the rest are zeros.
  }
  \item{G}{
	A square matrix representing the weighted graph structure of the predictors.
    Each entry \eqn{G[i, j]} indicates the weight of the edge between node \eqn{i} and node \eqn{j}.
    A value of zero means no edge is present. Diagonal entries are ignored.
  }
}
\details{
  The design covariance matrix \eqn{\Sigma} is first constructed according to the argument
	\code{xcov}, and the corresponding precision matrix is computed as
	\eqn{\Theta = \Sigma^{-1}}. The matrix \code{G} is then obtained by thresholding the
	entries of \eqn{\Theta}: \eqn{G[i, j]} is set to \eqn{|\,\Theta_{ij}\,|} if 
	\eqn{|\,\Theta_{ij}\,| > 10^{-8}} and to \eqn{0} otherwise.
	Thus, \code{G} can be interpreted as a weighted adjacency matrix; if only the graph
	structure is of interest, it can be further binarized by treating all nonzero entries
	as weight 1.
}
\references{
Tan, X., Zhang, X., Cui, Y., and Liu, X. (2024).
Uncertainty quantification in high-dimensional linear models incorporating graphical structures with applications to gene set analysis.
\emph{Bioinformatics}, \bold{40}(9), btae541.

Cai, L., Guo, X., Lian, H., and Zhu, L. (2025).
Statistical inference for high-dimensional convoluted rank regression.
\emph{Journal of the American Statistical Association}, 1--25.

Zhong, P.-S. and Chen, S. X. (2011).
Tests for high-dimensional regression coefficients with factorial designs.
\emph{Journal of the American Statistical Association}, \bold{106}(493), 260--274.
}
\examples{
set.seed(123)
dat <- Data_simu(n = 50, p = 100, s0 = 5,
                 xcov = "AR1", error = "norm", signal = 1)
str(dat)
} 
