\name{QRMissingUniBayesMix}
\alias{QRMissingUniBayesMix}
\alias{coef.QRMissingUniBayesMix}
\alias{plot.QRMissingUniBayesMix}
\alias{summary.QRMissingUniBayesMix}
\title{Quantile Regression in the Presence of Monotone Missingness; Univariate case with mixture normals; using HMC}
\usage{
QRMissingUniBayesMix(y, R, X, tau = 0.5, mcmc, prior)

\method{coef}{QRMissingUniBayesMix}(mod, ...)

\method{summary}{QRMissingUniBayesMix}(mod, ...)

\method{plot}{QRMissingUniBayesMix}(mod, ...)
}
\arguments{
  \item{y}{}

  \item{R}{}

  \item{X}{}

  \item{tau}{}

  \item{mcmc}{}

  \item{prior}{}
}
\value{
An object of class \code{QRMissingUni}, for which many
generic functions are available.
}
\description{
Function for quantile regression in the presence of
monotone missingness for univariate case
}
\details{
Using Bayesian method for univariate case for mixture of
normals
}
\author{
Minzhao Liu, Mike Daniels
}

