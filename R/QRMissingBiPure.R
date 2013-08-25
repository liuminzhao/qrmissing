##' Function for quantile regression in the presence of monotone
##' missingness for bivariate case using \emph{Pure R}
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
QRMissingBiPure <- function(y, R, X, tau = 0.5, sp = NULL,
                            init = NULL, method = 'uobyqa',
                            tol = 0.00001, control = list(maxit = 1000,
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

  ## negative log likelihood function
  nll <- function(param){
    ## translate param
    gamma1 <- param[1:xdim]
    beta1 <- param[(xdim + 1):(2*xdim)]
    sigma11 <- exp(param[3*xdim + 3])
    sigma10 <- exp(param[3*xdim + 4])
    gamma2 <- param[(2*xdim + 1):(3*xdim)]
    beta2 <- sp[1:xdim] # SP for R = 0
    sigma21 <- exp(param[3*xdim + 5])
    sigma20 <- exp(param[3*xdim + 5] + sp[xdim + 2])
    betay <- param[3*xdim + 1] # for R = 1
    beta2y <- betay + sp[xdim + 1] # SP for R = 0
    p <- exp(param[3*xdim + 2])/(1 + exp(param[3*xdim + 2]))

    ## Delta1 function
    Delta1 <- function(x){
      quan <- gamma1 %*% x
      lp <- beta1 %*% x
      sigma1 <- sigma11
      sigma0 <- sigma10
      target1 <- function(d){
        return(tau - p*pnorm((quan - d - lp)/sigma1) - (1 - p)*pnorm(
          (quan - d + lp)/sigma0))
      }
      interval <- c(-10, 10)

      ans <- uniroot.all(target1, interval, tol = tol)[1]
      rootiter <- 0
      repeat {
        if (!is.na(ans)) {
          break
        } else {
          interval <- interval * 2
          ans <- uniroot.all(target1, interval, tol = tol)[1]
        }
        rootiter <- rootiter + 1
        if (rootiter > 50) {
          cat('can not bracket root fot d1 \n')
          break
        }
      }
      return(ans)
    }

    Delta2 <- function(x){
      d1 <- Delta1(x)
      quan1 <- gamma1 %*% x
      lp1 <- beta1 %*% x
      sigma1 <- sigma11
      sigma0 <- sigma10
      mu11 <- d1 + lp1
      mu10 <- d1 - lp1
      quan2 <- gamma2 %*% x
      lp2 <- beta2 %*% x
      sigma2 <- sigma21
      sigma2sp <- sigma20

      target2 <- function(d2){
        if (betay == 0){
          int1 <- pnorm((quan2 - d2)/sigma2)
        } else {
          int1 <- pnorm(((quan2 - d2)/betay - mu11)/sqrt(sigma2^2/betay^2 + sigma1^2))
          if (betay < 0) {
            int1 <- 1 - int1
          }
        }

        if (beta2y == 0){
          int2 <- pnorm((quan2 - d2 - lp2)/sigma2sp)
        } else {
          int2 <- pnorm(((-d2 - lp2 + quan2)/beta2y - mu10)/sqrt(sigma2sp^2/beta2y^2 + sigma0^2))
          if (beta2y < 0){
            int2 <- 1 - int2
          }
        }
        return(tau - p*int1 - (1-p)*int2)
      }
      interval <- c(-10, 10)
      ans <- c(d1, uniroot.all(target2, interval, tol = tol)[1])
      rootiter <- 0
      repeat {
        if (!is.na(ans[2])) {
          break
        } else {
          interval <- interval * 2
          ans <- c(d1, uniroot.all(target2, interval, tol = tol)[1])
        }
        rootiter <- rootiter + 1
        if (rootiter > 50) {
          cat('can not bracket the root for d2 \n')
          break
        }
      }
      return(ans)
    }

    ## get delta
    delta <- t(apply(X, 1, function(l) Delta2(l)))
    d1 <- delta[,1]
    d2 <- delta[,2]

    ## Y1
    lp1 <- X %*% as.matrix(beta1)
    mu11 <- d1 + lp1
    mu10 <- d1 - lp1
    sigma1 <- sigma11
    sigma0 <- sigma10

    ## Y2
    mu21 <- d2 + betay * y[,1]
    sigma2 <- sigma21

    ## ll
    ll11 <- sum(dnorm(y[,1], mu11, sigma1, log=T)[R==1])
    ll10 <- sum(dnorm(y[,1], mu10, sigma0, log=T)[R==0])
    ll21 <- sum(dnorm(y[,2], mu21, sigma2, log=T)[R==1])
    ans <- ll11+ll10+ll21+ num*log(p) + (n - num)*log(1 - p)

    return(-ans)

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
##   mod$d <- d

  class(mod) <- "QRMissingBi"

  return(mod)

}
