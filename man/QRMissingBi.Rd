\name{QRMissingBi}
\alias{QRMissingBi}
\alias{coef.QRMissingBi}
\alias{plot.QRMissingBi}
\alias{print.QRMissingBi}
\alias{summary.QRMissingBi}
\title{Quantile Regression in the Presence of Monotone Missingness; Bivariates case}
\usage{
QRMissingBi(formula, R, tau = 0.5, sp = NULL, init = NULL,
  method = "uobyqa", control = list(maxit = 1000, trace = 0),
  hess = FALSE, model = "slope")

\method{coef}{QRMissingBi}(mod, ...)

\method{print}{QRMissingBi}(mod, ...)

\method{summary}{QRMissingBi}(mod, ...)

\method{plot}{QRMissingBi}(mod, ...)
}
\arguments{
  \item{y}{[n, 2] response.}

  \item{R}{[n]: 0-1 missingness indicator where \eqn{R = 1}
  indicates both \eqn{Y_1} and \eqn{Y_2} are observed,
  while \eqn{R = 0} indicates \eqn{Y_2} is missing}

  \item{X}{[n, p]: covariates matrix including intercept}

  \item{tau}{: quantile requested. Only single quantile is
  allowed for now.}

  \item{sp}{: [p + 2] sensitivity paramters. \code{sp[1:p]}
  is \eqn{\beta_2} and \code{sp[p+1]} is the sensitivity
  parameter for \eqn{\beta_y^(0)}, \code{sp[p+2]} is the
  one for \eqn{\sigma_2^(0)}}

  \item{init}{: initial value for all parameters. The
  default are estimates by \code{rq} function from
  \pkg{quantreg} and sample proportion of missingness.}

  \item{method}{: optimization method, default: 'uobyqa'
  from minqa. Other allowed methods are 'BFGS',
  'Nelder-Mead', BFGS', 'CG', 'L-BFGS-B', 'Nelder-Mead')}

  \item{control}{: optimzation control option}

  \item{hess}{: whether to calculate hessian matrix,
  default is FALSE}
}
\value{
An object of class \code{QRMissingBi}, for which many
generic functions are available.
}
\description{
Function for quantile regression in the presence of
monotone missingness for bivariate case

Plot the fitted residuals to check the normal assumption
for goodness of fit
}
\details{
Fit a quantile regression for bivariate responses. For more
details, please check the paper.
}
\author{
Minzhao Liu, Mike Daniels
}

