\name{infer_norm}
\alias{infer_norm}
\title{
    Normal-approximation inference for graph-assisted convoluted rank regression
}
\description{
    Perform coordinatewise inference for high-dimensional linear models fitted
    via graph-assisted convoluted rank regression. Given an initial estimator,
    its debiased version, and a graph structure \code{G}, the function constructs
    marginal confidence intervals for each coefficient using a normal approximation.
    It also reports standard error estimates and two-sided p-values. When the true
    coefficients are provided, the function additionally evaluates empirical coverage
    in simulation studies.
}
\usage{
    infer_norm(X, y, beta_hat, beta_true = NULL, debias, G, h, alpha = 0.05, 
      inv_method = "Tan2024")
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
  \item{inv_method}{
  A character string specifying the method used to estimate the inverse Gram
  (or precision) matrix in the variance formula. Possible values are:
  \describe{
    \item{\code{"Tan2024"}}{
      The graph-assisted inverse Gram matrix estimator proposed by
      Tan, Zhang, Cui, and Liu (2024). In the implementation, the inverse Gram
      matrix is computed via the function
      \code{calculate_inv_gram(X, G)} from the \pkg{GCDL} package, which
      incorporates the graph structure \code{G}.
    }
    \item{\code{"Cai2025"}}{
      The precision matrix estimator used in Cai, Guo, Lian, and Zhu (2025).
      In the implementation, a sparse precision matrix is estimated using the
      function \code{sugm} from the \pkg{flare} package.
    }
  }
  Default is \code{"Tan2024"}.
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
  \item{Inv_Gram}{
    The estimated inverse Gram/precision matrix used in the inference procedure.
  }
  \item{sehat}{
    A numeric vector of length \eqn{p} giving the estimated standard errors for
    each coordinate of the debiased estimator.
  }
  \item{p_value}{
  A numeric vector of length \eqn{p} giving two-sided p-values based on the
  normal approximation for testing the null hypotheses
  \eqn{H_0: \beta_j = 0}, \eqn{j = 1, \dots, p}.
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

## Normal-approximation inference
res <- infer_norm(X, y, beta_hat,
                  beta_true = beta_true,
                  debias = beta_db,
                  G = G, h = 1,
                  inv_method = "Tan2024")
}
