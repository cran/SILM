% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Sim_CI.R
\name{Sim.CI}
\alias{Sim.CI}
\title{Simultaneous Confidence Interval}
\usage{
Sim.CI(X, Y, set, M = 500, alpha = 0.95)
}
\arguments{
\item{X}{n times p design matrix.}

\item{Y}{Response variable.}

\item{set}{The set of variables of interest.}

\item{M}{The number of bootstrap replications (default 500).}

\item{alpha}{The nominal level alpha (default 0.95).}
}
\value{
The de-biased Lasso estimator, the confidence bands (lower bound and
 upper bound) delivered by the non-studentized and the studentized statistics.
}
\description{
This function implements the method for constructing simultaneous confidence
interval in Zhang and Cheng (2017).
}
\references{
Zhang, X., and Cheng, G. (2017) Simultaneous Inference for High-dimensional
Linear Models, \emph{Journal of the American Statistical Association}, 112, 757-768.
}
\examples{
## The function is intended for large n and p.
## Use small p here for illustration purpose only.
n <- 100
p <- 10
s0 <- 3
set <- 1:s0
Sigma <- matrix(NA, p, p)
for (i in 1:p) Sigma[i,] <- 0.9^(abs(i-(1:p)))
X <- matrix(rnorm(n*p), n, p)
X <- t(t(chol(Sigma))\%*\%t(X))
beta <- rep(0,p)
beta[1:s0] <- runif(s0,0,2)
Y <- X\%*\%beta+rt(n,4)/sqrt(2)
Sim.CI(X, Y, set)
}
