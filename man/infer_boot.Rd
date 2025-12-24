\name{infer_boot}
\alias{infer_boot}
\title{
    Bootstrap-based inference for graph-assisted convoluted rank regression
}
\description{
    Perform coordinatewise inference for high-dimensional linear models fitted
    via graph-assisted convoluted rank regression. Given an initial estimator,
    its debiased version, and a graph structure \code{G}, the function applies
    a multiplier bootstrap based on the smoothed rank loss to construct
    marginal confidence intervals for each coefficient. When the true coefficients
    are provided, the function additionally evaluates empirical coverage in
    simulation studies.
}
\usage{
    infer_boot(X, y, beta_hat, beta_true = NULL, debias, G, h, alpha = 0.05, B = 1000)
}
\arguments{
  \item{X}{
    An \eqn{n \times p} design matrix of predictors.
  }
  \item{y}{
    A numeric response vector of length \eqn{n}.
  }
  \item{beta_hat}{
    A numeric vector of length \eqn{p} giving the initial estimator of the
    regression coefficients (e.g., the output from \code{Lfabs} or
    \code{ebic.Lfabs}).
  }
  \item{beta_true}{
    Optional. A numeric vector of length \eqn{p} containing the true regression
    coefficients. If provided, it is used to assess empirical coverage of the
    confidence intervals in simulation studies; if \code{NULL} (default), the
    coverage component is not returned.
  }
  \item{debias}{
    A numeric vector of length \eqn{p} giving the debiased estimator of the
    regression coefficients, for example from \code{debias_lap}.
  }
  \item{G}{
	  A square matrix representing the weighted graph structure of the predictors.
    Each entry \eqn{G[i, j]} indicates the weight of the edge between node \eqn{i} and node \eqn{j}.
    A value of zero means no edge is present. Diagonal entries are ignored.
  }
  \item{h}{
    A positive bandwidth parameter used in the smoothed (convoluted) rank loss
    and its derivatives. Default is \code{1}. It should match the bandwidth used in estimation.
  }
  \item{alpha}{
    Nominal significance level for the marginal confidence intervals.
    Default is \code{0.05}.
  }
  \item{B}{
    Number of bootstrap replications used in the algorithm.
    Default is \code{1000}.
  }
}
\value{
A list with the following components:
  \item{cover}{
    If \code{beta_true} is provided, a logical vector of length \eqn{p}. The \eqn{j}-th
    entry is \code{TRUE} if the true coefficient \eqn{\beta^{\ast}_j} is contained
    in the constructed confidence interval, and \code{FALSE} otherwise. If
    \code{beta_true = NULL}, this component is \code{NULL}.
  }
  \item{length}{
    A numeric vector of length \eqn{p} giving the lengths of the marginal
    confidence intervals for each coefficient.
  }
  \item{lower}{
    A numeric vector of length \eqn{p} giving the lower endpoints of the
    marginal confidence intervals.
  }
  \item{upper}{
    A numeric vector of length \eqn{p} giving the upper endpoints of the
    marginal confidence intervals.
  }
}
\references{
Tan, X., Zhang, X., Cui, Y., and Liu, X. (2024).
Uncertainty quantification in high-dimensional linear models incorporating graphical structures with applications to gene set analysis.
\emph{Bioinformatics}, \bold{40}(9), btae541.

Cai, L., Guo, X., Lian, H., and Zhu, L. (2025).
Statistical inference for high-dimensional convoluted rank regression.
\emph{Journal of the American Statistical Association}, 1--25.
}

\examples{
set.seed(123)

## Simulate data (see Data_simu for details)
dat <- Data_simu(n = 50, p = 100, s0 = 5,
                 xcov = "AR1", error = "norm", signal = 1)
X <- dat$X
y <- dat$y
G <- dat$G
beta_true <- dat$beta_T

## Initial estimator
fit_init <- ebic.Lfabs(X, y, G, h = 1)
beta_hat <- fit_init$beta

## Debiased estimator
beta_db <- debias_lap(X, y, beta_hat, G, h = 1)

## Bootstrap-based inference
res <- infer_boot(X, y, beta_hat, beta_true, beta_db, G, h = 1)
}
