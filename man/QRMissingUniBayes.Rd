\name{QRMissingUniBayes}
\alias{coef.QRMissingUniBayes}
\alias{plot.QRMissingUniBayes}
\alias{QRMissingUniBayes}
\alias{summary.QRMissingUniBayes}
\title{Quantile Regression in the Presence of Monotone Missingness; Univariate case}
\usage{
  QRMissingUniBayes(y, R, X, tau = 0.5, mcmc, prior)

  \method{coef}{QRMissingUniBayes} (mod, ...)

  \method{summary}{QRMissingUniBayes} (mod, ...)

  \method{plot}{QRMissingUniBayes} (mod, ...)
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
  Using Bayesian method for univariate case for single
  normal error
}
\author{
  Minzhao Liu, Mike Daniels
}

