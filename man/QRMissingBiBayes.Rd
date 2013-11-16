\name{QRMissingBiBayes}
\alias{coef.QRMissingBiBayes}
\alias{plot.QRMissingBiBayes}
\alias{QRMissingBiBayes}
\alias{summary.QRMissingBiBayes}
\title{Quantile Regression in the Presence of Monotone Missingness; Bivariate case}
\usage{
  QRMissingBiBayes(y, R, X, tau = 0.5, mcmc, prior)

  \method{coef}{QRMissingBiBayes} (mod, ...)

  \method{summary}{QRMissingBiBayes} (mod, ...)

  \method{plot}{QRMissingBiBayes} (mod, ...)
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
  An object of class \code{QRMissingBiBayes}, for which
  many generic functions are available.
}
\description{
  Function for quantile regression in the presence of
  monotone missingness for Bivariate case
}
\details{
  Using Bayesian method for Bivariate quantile regression
  with single normal error assumption
}
\author{
  Minzhao Liu, Mike Daniels
}

