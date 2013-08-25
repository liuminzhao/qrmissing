\name{residuals}
\alias{residuals}
\title{residuals for fitted value}
\usage{
  residuals(param, y, X, R, tau, sp)
}
\arguments{
  \item{param}{}

  \item{y}{responses}

  \item{X}{covariates matrix}

  \item{R}{missingnes indicator}

  \item{tau}{quantile of interest}

  \item{sp}{sensitivity parameters}
}
\value{
  a n by 2 matrix of residuals
}
\description{
  Residuals for fitted value
}
\details{
  Residuals for fitted value after pluging MLE
}
\author{
  Minzhao Liu
}

