\name{QRMissingBiBayesMix}
\alias{QRMissingBiBayesMix}
\alias{coef.QRMissingBiBayesMix}
\alias{plot.QRMissingBiBayesMix}
\alias{summary.QRMissingBiBayesMix}
\title{Quantile Regression in the Presence of Monotone Missingness; Bivariate case}
\usage{
QRMissingBiBayesMix(y, R, X, tau = 0.5, mcmc, prior)

\method{coef}{QRMissingBiBayesMix}(mod, ...)

\method{summary}{QRMissingBiBayesMix}(mod, ...)

\method{plot}{QRMissingBiBayesMix}(mod, ...)
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
An object of class \code{QRMissingBiBayes}, for which many
generic functions are available.
}
\description{
Function for quantile regression in the presence of
monotone missingness for Bivariate case
}
\details{
Using Bayesian method for Bivariate quantile regression
with mixture normal error assumption
}
\author{
Minzhao Liu, Mike Daniels
}

