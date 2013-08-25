##' Function for quantile regression in the presence of monotone
##' missingness for bivariate case
##'
##' Fit a quantile regression for bivariate responses. For more
##' details, please check the paper.
##'
##' @title Quantile Regression in the Presence of Monotone Missingness; Bivariates case
##' @param y [n, 2] response.
##' @param R [n]: 0-1 missingness indicator where \eqn{R = 1}
##' indicates both \eqn{Y_1} and \eqn{Y_2} are observed, while \eqn{R
##' = 0} indicates \eqn{Y_2} is missing
##' @param X [n, p]: covariates matrix including intercept
##' @param tau : quantile requested. Only single quantile is allowed
##' for now.
##' @param sp : [p + 2] sensitivity paramters. \code{sp[1:p]} is
##' \eqn{\beta_2} and \code{sp[p+1]} is the sensitivity parameter for
##' \eqn{\beta_y^(0)}, \code{sp[p+2]} is the one for \eqn{\sigma_2^(0)}
##' @param init : initial value for all parameters. The default are
##' estimates by \code{rq} function from \pkg{quantreg} and sample
##' proportion of missingness.
##' @param method : optimization method, default: 'uobyqa' from minqa. Other allowed methods are 'BFGS', 'Nelder-Mead', BFGS', 'CG', 'L-BFGS-B', 'Nelder-Mead')
##' @param control : optimzation control option
##' @param hess : whether to calculate hessian matrix, default is FALSE
##' @return An object of class \code{QRMissingBi}, for which many
##' generic functions are available.
##' @author Minzhao Liu
##' @export
QRMissingBi <- function(y, R, X, tau = 0.5, sp = NULL,
                        init = NULL, method = 'uobyqa',
                        control = list(maxit = 1000,
                          trace = 0), hess = FALSE){
  ## data
  n <- dim(y)[1]
  num <- sum(R)
  xdim <- dim(X)[2]

  ## initial
  if (is.null(sp)) {
    sp <- rep(0, xdim  + 2)
  }
  if (!is.null(init)){
    param <- init
  } else {
    lmcoef1 <- coef(rq(y[,1] ~ X[,-1], tau = tau))
    lmcoef2 <- coef(rq(y[,2][R == 1] ~ X[R == 1,-1], tau = tau))
    param <- rep(0, 3*xdim + 5)
    param[1:xdim] <- lmcoef1
    param[(2*xdim + 1):(3*xdim)] <- lmcoef2
    param[3*xdim + 2] = log(num/(n-num))
  }

  ## nll
  nll <- function(param){
    ll2(param, y, X, R, tau, sp)
  }

  ## optimize nll to get MLE
  optim_method <- c('BFGS', 'CG', 'L-BFGS-B', 'Nelder-Mead')

  if (method %in% optim_method) {
    mod <- optim(param, nll, method = method, control = control)
  } else {
    minqa_control <- list(iprint = control$trace, maxfun = 20*control$maxit)
    if (method == 'bobyqa'){
      mod <- bobyqa(param, nll, control=minqa_control)
    } else if (method == 'uobyqa') {
      mod <- uobyqa(param, nll, control=minqa_control)
    } else if (method == 'newuoa') {
      mod <- newuoa(param, nll, control=minqa_control)
    }
  }

  ## residuals
  res <- residuals(mod$par, y, X, R, tau, sp)

  ## Hessian matrix and grad
  if (hess) {
    Hessian <- hessian(nll, mod$par)
    d <- grad(nll, mod$par)
    Jninv <- solve(Hessian)
    se <- matrix(0, 2, xdim)
    se[1, ] <- sqrt(diag(Jninv)[1:xdim])
    se[2, ] <- sqrt(diag(Jninv)[(2*xdim + 1):(3*xdim)])
    rownames(se) <- c('Q1', 'Q2')
  } else {
    Hessian <- NULL
    se <- NULL
  }

  mod$n <- n
  mod$xdim <- xdim
  mod$X <- X
  mod$y <- y
  mod$R <- R
  mod$tau <- tau
  mod$method <- method
  mod$Hessian <- Hessian
  mod$se <- se
  mod$res <- res

  class(mod) <- "QRMissingBi"

  return(mod)

}

##' @rdname QRMissingBi
##' @method coef QRMissingBi
##' @S3method coef QRMissingBi
coef.QRMissingBi <- function(mod, ...){
  q <- mod$xdim
  param <- mod$par[c(1:q, (2*q + 1):(3*q))]
  coef <- matrix(param, 2, q, byrow = T)
  rownames(coef) <- c('Q1Coef', 'Q2Coef')
  return(coef)
}

##' @rdname QRMissingBi
##' @method print QRMissingBi
##' @S3method print QRMissingBi
print.QRMissingBi <- function(mod, ...){
  cat('Coefficients: \n')
  print(coef(mod))
}

##' @rdname QRMissingBi
##' @method summary QRMissingBi
##' @S3method summary QRMissingBi
summary.QRMissingBi <- function(mod, ...){
  n <- mod$n
  R <- mod$R
  tau <- mod$tau
  param <- mod$par
  q <- mod$xdim

  cat('Number of observations: ', n, '\n')
  cat('Sample proportion of observed data: ', sum(R)/n, '\n')
  cat('Estimated pi:', exp(param[3*q + 2])/(1 + exp(param[3*q + 2])), '\n')
  cat('Quantile: ', tau, '\n')
  cat('Optimization method: ', mod$method, '\n')
  optim_method <- c('BFGS', 'CG', 'L-BFGS-B', 'Nelder-Mead')

  if (mod$method %in% optim_method) {
    cat('Model converged: ', ifelse(mod$convergence, 'No', 'Yes'), '\n')
  } else {
    cat('Model converged: ', ifelse(mod$ierr, 'No', 'Yes'), '\n')
  }
  cat('Quantile regression coefficients: \n')
  print(coef(mod))
  cat('Standard error: \n')
  print(mod$se)

}

##' A goodness of fit check
##'
##' Plot the fitted residuals to check the normal assumption for
##' goodness of fit
##'
##' @rdname QRMissingBi
##' @method plot QRMissingBi
##' @S3method plot QRMissingBi
plot.QRMissingBi <- function(mod, ...){
  R <- mod$R
  res <- mod$res
  tau <- mod$tau
  par(mfrow = c(1, 2))
  qqnorm(res[, 1], main = bquote(paste('Normal Q-Q Plot for ', Y[1], ' of quantile ', tau == .(tau))))
  qqline(res[, 1])
  qqnorm(res[R == 1, 1], main = bquote(paste('Normal Q-Q Plot for ', Y[2], ' of quantile ', tau == .(tau))))
  qqline(res[R == 1, 1])
}

##' Observed Negative Log Likelihood
##'
##' Give the observed nagetive log likelihood
##' @title observed negative log likelihood
##' @param param
##' @param y
##' @param X
##' @param R
##' @param tau
##' @param sp
##' @return negative log likelihood
##' @useDynLib qrmissing
##' @author Minzhao Liu
ll2 <- function(param, y, X, R, tau, sp){
  n <- dim(y)[1]
  xdim <- dim(X)[2]
  num <- sum(R)

  gamma1 <- param[1:xdim]
  beta1 <- param[(xdim + 1):(2*xdim)]
  sigma1 <- c(exp(param[3*xdim + 3]), exp(param[3*xdim + 4]))
  gamma2 <- param[(2*xdim + 1):(3*xdim)]
  beta2sp <- sp[1:xdim] # SP for R = 0
  sigma21 <- exp(param[3*xdim + 5])
  sigma21sp <- exp(sp[xdim + 2])  # SP for R = 0
  betay <- param[3*xdim + 1] # for R = 1
  betaysp <- sp[xdim + 1] # SP for R = 0
  p <- exp(param[3*xdim + 2])/(1 + exp(param[3*xdim + 2]))

  d <- matrix(0, n, 2)
  d <- .Fortran("mydelta2",
                x = as.double(X),
                gamma1 = as.double(gamma1),
                beta1 = as.double(beta1),
                sigma1 = as.double(sigma1),
                gamma2 = as.double(gamma2),
                beta2sp = as.double(beta2sp),
                sigma21 = as.double(sigma21),
                sigma21sp = as.double(sigma21sp),
                betay = as.double(betay),
                betaysp = as.double(betaysp),
                p = as.double(p),
                tau = as.double(tau),
                n = as.integer(n),
                xdim = as.integer(xdim),
                delta = as.double(d))$delta

  d <- matrix(d, n, 2)

  lp1 <- X %*% beta1
  mu11 <- d[, 1] + lp1
  mu10 <- d[, 1] - lp1
  mu21 <- d[, 2] + betay * y[, 1]
  ll11 <- sum(dnorm(y[, 1], mu11, sigma1[1], log=T)[R==1])
  ll10 <- sum(dnorm(y[, 1], mu10, sigma1[2], log=T)[R==0])
  ll21 <- sum(dnorm(y[, 2], mu21, sigma21, log = T)[R==1])
  ans <- ll11 + ll10 + ll21 + num*log(p) + (n - num)*log(1 - p)

  return(-ans)
}

##' Residuals for fitted value
##'
##' Residuals for fitted value after pluging MLE
##'
##' @title residuals for fitted value
##' @param param
##' @param y responses
##' @param X covariates matrix
##' @param R missingnes indicator
##' @param tau quantile of interest
##' @param sp sensitivity parameters
##' @return a n by 2 matrix of residuals
##' @author Minzhao Liu
residuals <- function(param, y, X, R, tau, sp){
  n <- dim(y)[1]
  xdim <- dim(X)[2]
  num <- sum(R)

  gamma1 <- param[1:xdim]
  beta1 <- param[(xdim + 1):(2*xdim)]
  sigma1 <- c(exp(param[3*xdim + 3]), exp(param[3*xdim + 4]))
  gamma2 <- param[(2*xdim + 1):(3*xdim)]
  beta2sp <- sp[1:xdim] # SP for R = 0
  sigma21 <- exp(param[3*xdim + 5])
  sigma21sp <- exp(sp[xdim + 2])  # SP for R = 0
  betay <- param[3*xdim + 1] # for R = 1
  betaysp <- sp[xdim + 1] # SP for R = 0
  p <- exp(param[3*xdim + 2])/(1 + exp(param[3*xdim + 2]))

  d <- matrix(0, n, 2)
  d <- .Fortran("mydelta2",
                x = as.double(X),
                gamma1 = as.double(gamma1),
                beta1 = as.double(beta1),
                sigma1 = as.double(sigma1),
                gamma2 = as.double(gamma2),
                beta2sp = as.double(beta2sp),
                sigma21 = as.double(sigma21),
                sigma21sp = as.double(sigma21sp),
                betay = as.double(betay),
                betaysp = as.double(betaysp),
                p = as.double(p),
                tau = as.double(tau),
                n = as.integer(n),
                xdim = as.integer(xdim),
                delta = as.double(d))$delta

  d <- matrix(d, n, 2)

  lp1 <- X %*% beta1
  mu11 <- d[, 1] + lp1
  mu10 <- d[, 1] - lp1
  mu21 <- d[, 2] + betay * y[, 1]

  res <- matrix(NA, n, 2)
  res[R == 1, 1] <- (y[R==1, 1] - mu11[R == 1])/sigma1[1]
  res[R == 0, 1] <- (y[R==0, 1] - mu10[R == 0])/sigma1[2]
  res[R == 1, 2] <- (y[R==1, 2] - mu21[R == 1])/sigma21

  return(res)
}
