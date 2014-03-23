##' Function for quantile regression in the presence of monotone
##' missingness for univariate case
##'
##' Using Bayesian method for univariate case
##' for single normal error
##'
##' @title Quantile Regression in the Presence of Monotone Missingness; Univariate case
##' @param y
##' @param R
##' @param X
##' @param tau
##' @param mcmc
##' @param prior
##' @return An object of class \code{QRMissingUni}, for which many
##' generic functions are available.
##' @author Minzhao Liu, Mike Daniels
##' @export
QRMissingUniBayes <- function(y, R, X, tau = 0.5,
                              mcmc, prior
                              ){

  ## data
  n <- length(y)
  xdim <- dim(X)[2]
  num <- sum(R)

  ## prior
  betapm <- prior$betapm
  betapv <- prior$betapv
  gammapm <- prior$gammapm
  gammapv <- prior$gammapv
  a <- prior$a ## for sigma
  b <- prior$b
  c <- prior$c ## for p
  d <- prior$d

  ## MCMC
  nsave <- mcmc$nsave
  nskip <- mcmc$nskip
  nburn <- mcmc$nburn
  ndisp <- mcmc$ndisp

  ## SAVE
  betasave <- gammasave <- matrix(0, nsave, xdim)
  psave <- rep(0, nsave)
  sigmasave <- matrix(0, nsave, 2)

  ## TUNE
  tunegamma <- tunebeta <- 0.3
  tunep <- 0.1
  tunesigma <- 0.3
  arate <- 0.25
  att <- acc <- 0

  ## initial

  gamma <- beta <- rep(0, xdim)
  p <- num/n
  sigma <- rep(1, 2)

  isave <- 0
  skipcount <- 0
  dispcount <- 0
  nscan <- nburn + nskip * nsave

  start <- proc.time()[3]

########################################

  ## first

  loglikeo <- LLUniSingle(gamma, beta, p, sigma, tau, y, X, R)

  ## roll

  for (iscan in 1:nscan) {

    ## GAMMA
    att <- att + 1
    gammac <- rnorm(xdim, gamma, tunegamma)

    ## BETA
    betac <- rnorm(xdim, beta, tunebeta)

    ## P
    pc <- rnorm(1, p, tunep)
    pc <- max(min(pc, 0.99), 0.01)

    ## SIGMA
    theta <- log(sigma)
    thetac <- rnorm(2, theta, tunesigma)
    sigmac <- exp(thetac)

    ## LOG LIKELIHOOD CANDIDATE
    loglikec <- LLUniSingle(gammac, betac, pc, sigmac, tau, y, X, R)

    ## PRIOR
    logpriorc <- sum(dnorm(gammac, gammapm, gammapv, log = T))
    logprioro <- sum(dnorm(gamma, gammapm, gammapv, log = T))
    logpriorc <- logpriorc + sum(dnorm(betac, betapm, betapv, log = T))
    logprioro <- logprioro + sum(dnorm(beta, betapm, betapv, log = T))
    logpriorc <- logpriorc + dbeta(pc, c/2, d/2, log = T)
    logprioro <- logprioro + dbeta(p, c/2, d/2, log = T)
    logpriorc <- logpriorc + sum(dgamma(sigmac, a/2, b/2, log = T))
    logprioro <- logprioro + sum(dgamma(sigma, a/2, b/2, log = T))
    logcgkc <- -sum(theta)
    logcgko <- -sum(thetac)

    ## ACCEPTA

    ratio <- loglikec + logpriorc - loglikeo - logprioro + logcgkc - logcgko

    if (log(runif(1)) <= ratio) {
      acc <- acc + 1
      loglikeo <- loglikec
      gamma <- gammac
      beta <- betac
      p <- pc
      sigma <- sigmac
    }

    ## TUNE
    if (att >= 100 && iscan < nburn) {
      prop <- acc/att
      tunegamma <- tunegamma*ifelse(prop > arate, 2, 0.5)
      tunebeta <- tunebeta*ifelse(prop > arate, 2, 0.5)
      tunesigma <- tunesigma*ifelse(prop > arate, 2, 0.5)
      tunep <- tunep*ifelse(prop > arate, 2, 0.5)
      att <- acc <- 0
    }

    ## SAVE

    if (iscan > nburn) {
      skipcount <- skipcount + 1
      if (skipcount >= nskip) {
        isave <- isave + 1
        dispcount <- dispcount + 1
        gammasave[isave, ] <- gamma
        betasave[isave, ] <- beta
        psave[isave] <- p
        sigmasave[isave, ] <- sigma
        skipcount <- 0
        if (dispcount >= ndisp) {
          dispcount <- 0
          cat(isave, 'out of', nsave, proc.time()[3] - start, '\n')
        }
      }
    }
  }

  ans <- list(gammasave = gammasave,
              betasave = betasave,
              psave = psave,
              sigmasave = sigmasave,
              n = n,
              xdim = xdim,
              y = y,
              X = X,
              R = R,
              tau = tau,
              tune = list(gamma = tunegamma, beta = tunebeta,
                sigma = tunesigma, p = tunep)
              )

  class(ans) <- 'QRMissingUniBayes'

  return(ans)

}

##' Observed Log Likelihood
##'
##'
##' @title observed log likelihood for univariate case with single normal
##' @param gamma
##' @param beta
##' @param p
##' @param sigma
##' @param tau
##' @param y
##' @param X
##' @param R
##' @return log likelihood of univariate response with missingness
##' @useDynLib qrmissing
##' @author Minzhao Liu, Mike Daniels
##' @export
LLUniSingle <- function(gamma, beta, p, sigma, tau, y, X, R){
  n <- length(y)
  xdim <- dim(X)[2]
  num <- sum(R)

  d <- rep(0, n)
  d <- .Fortran("mydelta1bise",
                x = as.double(X),
                gamma = as.double(gamma),
                beta = as.double(beta),
                sigma = as.double(sigma),
                p = as.double(p),
                tau = as.double(tau),
                n = as.integer(n),
                xdim = as.integer(xdim),
                delta = as.double(d))$delta

  lp <- X %*% beta
  mu11 <- d + lp
  mu10 <- d - lp
  ll11 <- sum(dnorm(y, mu11, sigma[1], log=T)[R==1])
  ll10 <- sum(dnorm(y, mu10, sigma[2], log=T)[R==0])
  ans <- ll11 + ll10 + num*log(p) + (n - num)*log(1 - p)

  return(ans)
}


##' @rdname QRMissingUniBayes
##' @method coef QRMissingUniBayes
##' @S3method coef QRMissingUniBayes
coef.QRMissingUniBayes <- function(mod, ...){
  gamma <- apply(mod$gammasave, 2, mean)
  beta <- apply(mod$betasave, 2, mean)
  sigma <- apply(mod$sigmasave, 2, mean)
  p <- mean(mod$psave)
  return(list(gamma = gamma, beta = beta, sigma = sigma, p = p))
}

##' @rdname QRMissingUniBayes
##' @method summary QRMissingUniBayes
##' @S3method summary QRMissingUniBayes
summary.QRMissingUniBayes <- function(mod, ...){
  n <- mod$n
  R <- mod$R
  tau <- mod$tau

  cat('Number of observations: ', n, '\n')
  cat('Sample proportion of observed data: ', sum(R)/n, '\n')
  cat('Estimated pi:', coef(mod)$p, '\n')
  cat('Quantile: ', tau, '\n')
  cat('Quantile regression coefficients: \n')
  print(coef(mod))
}

##' @rdname QRMissingUniBayes
##' @method plot QRMissingUniBayes
##' @S3method plot QRMissingUniBayes
plot.QRMissingUniBayes <- function(mod, ...){
  xdim <- mod$xdim
  for (i in 1:xdim){
    plot(ts(mod$gammasave[, i]), main = paste('gamma', i, sep = ''))
  }
  for (i in 1:xdim){
    plot(ts(mod$betasave[, i]), main = paste('beta', i, sep = ''))
  }
  plot(ts(mod$psave), main = 'p')
  for (i in 1:2){
    plot(ts(mod$sigmasave[, i]), main = paste('sigma', i, sep = ''))
  }
}
