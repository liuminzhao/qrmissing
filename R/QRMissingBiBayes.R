##' Function for quantile regression in the presence of monotone
##' missingness for Bivariate case
##'
##' Using Bayesian method for Bivariate quantile regression
##' with single normal error assumption
##'
##' @title Quantile Regression in the Presence of Monotone Missingness; Bivariate case
##' @param y
##' @param R
##' @param X
##' @param tau
##' @param mcmc
##' @param prior
##' @return An object of class \code{QRMissingBiBayes}, for which many
##' generic functions are available.
##' @author Minzhao Liu, Mike Daniels
##' @export
QRMissingBiBayes <- function(y, R, X, tau = 0.5,
                              mcmc, prior
                              ){

  ## data
  n <- dim(y)[1]
  xdim <- dim(X)[2]
  num <- sum(R)

  ## prior
  betapm <- prior$betapm
  betapv <- prior$betapv
  gammapm <- prior$gammapm
  gammapv <- prior$gammapv
  sigmaa <- prior$sigmaa ## for sigma
  sigmab <- prior$sigmab
  alpha1 <- prior$alpha1 ## for p
  alpha2 <- prior$alpha2

  ## SP prior
  beta2pm <- prior$beta2pm
  beta2pv <- prior$beta2pv
  sigma21sppm <- prior$sigma21sppm
  sigma21sppv <- prior$sigma21sppv
  betaysppm <- prior$betaysppm
  betaysppv <- prior$betaysppv

  ## MCMC
  nsave <- mcmc$nsave
  nskip <- mcmc$nskip
  nburn <- mcmc$nburn
  ndisp <- mcmc$ndisp

  ## SAVE
  beta1save <- rep(0, nsave)
  gamma1save <- gamma2save <- matrix(0, nsave, xdim)
  sigma1save <- rep(0, nsave)
  sigma21save <- rep(0, nsave)
  betaysave <- rep(0, nsave)
  psave <- rep(0, nsave)

  ## TUNE
  tunegamma1 <- prior$tunegamma1
  tunegamma2 <- prior$tunegamma2
  tunebeta1 <- prior$tunebeta1
  tunesigma1 <- tunesigma21 <- prior$tunesigma1
  tunebetay <- prior$tunebetay
  tunep <- prior$tunep
  arate <- prior$arate

  att <- acc <- 0

  ## initial
  gamma1 <- coef(rq(y[, 1] ~ X[, -1], tau = tau))
  beta1 <- 0
  gamma2 <- coef(rq(y[R==1,2] ~ X[R==1, -1], tau = tau))
  beta2sp <- 0 # SP
  sigma1 <- 1
  sigma21 <- 1
  sigma21sp <- 0 # SP
  betay <- betaysp <- 0
  p <- num/n

  isave <- 0
  skipcount <- 0
  dispcount <- 0
  nscan <- nburn + nskip * nsave

  start <- proc.time()[3]

########################################

  ## first

  loglikeo <- LLBiSingle(gamma1, beta1, sigma1, gamma2, beta2sp, sigma21, sigma21sp, betay, betaysp, p, tau, y, X, R)

  ## roll

  for (iscan in 1:nscan) {

      att <- att + 1
      gamma1c <- rnorm(xdim, gamma1, tunegamma1)
      beta1c <- rnorm(1, beta1, tunebeta1)
      theta1 <- log(sigma1)
      theta1c <- rnorm(1, theta1, tunesigma1)
      sigma1c <- exp(theta1c)
      gamma2c <- rnorm(xdim, gamma2, tunegamma2)
      theta2 <- log(sigma21)
      theta2c <- rnorm(1, theta2, tunesigma21)
      sigma21c <- exp(theta2c)
      betayc <- rnorm(1, betay, tunebetay)
      pc <- rnorm(1, p, tunep)
      pc <- max(min(pc, 0.99), 0.01)

      ## Prior
      logprioro <- logpriorc <- 0

      logpriorc <- sum(dnorm(gamma1c, gammapm, gammapv, log = T)) + dnorm(beta1c, betapm, betapv, log = T) + dgamma(sigma1c, sigmaa, sigmab, log = T) + sum(dnorm(gamma2c, gammapm, gammapv, log = T)) + dgamma(sigma21c, sigmaa, sigmab, log = T) + dnorm(betayc, betapm, betapv, log = T) + dbeta(pc, alpha1, alpha2, log = T)

      logprioro <- sum(dnorm(gamma1, gammapm, gammapv, log = T)) + dnorm(beta1, betapm, betapv, log = T) + dgamma(sigma1, sigmaa, sigmab, log = T) + sum(dnorm(gamma2, gammapm, gammapv, log = T)) + dgamma(sigma21, sigmaa, sigmab, log = T) + dnorm(betay, betapm, betapv, log = T) + dbeta(p, alpha1, alpha2, log = T)

      ## additional prior
      logcgkc <- -theta1 - theta2
      logcgko <- -theta1c - theta2c


      ## loglike

      loglikec <- LLBiSingle(gamma1c, beta1c, sigma1c, gamma2c, beta2sp, sigma21c, sigma21sp, betayc, betaysp, pc, tau, y, X, R)

      ratio <- loglikec + logpriorc - loglikeo - logprioro + logcgkc - logcgko

      if (log(runif(1)) <= ratio) {
          acc <- acc + 1
          loglikeo <- loglikec
          gamma1 <- gamma1c
          beta1 <- beta1c
          sigma1 <- sigma1c
          gamma2 <- gamma2c
          sigma21 <- sigma21c
          betay <- betayc
          p <- pc
      }

      ## SP
      beta2sp <- rnorm(1, beta2pm, beta2pv)
      sigma21sp <- 0
      betaysp <- 0

    ## TUNE
      ## if (att >= 100 & iscan < nburn) {
      if (FALSE) {
          tunegamma1 <- pmin(pmax(tunegamma1*ifelse(acc/att > arate, 1.5, 0.7), 0.01), 10)
          tunebeta1 <- pmin(pmax(tunebeta1*ifelse(acc/att > arate, 1.5, 0.7), 0.01), 10)
          tunesigma1 <- pmin(pmax(tunesigma1*ifelse(acc/att > arate, 1.5, 0.7), 0.01), 10)
          tunegamma2 <- pmin(pmax(tunegamma2*ifelse(acc/att > arate, 1.5, 0.7), 0.01), 10)
          tunesigma21 <- pmin(pmax(tunesigma21*ifelse(acc/att > arate, 1.5, 0.7), 0.01), 10)
          tunebetay <- pmin(pmax(tunebetay*ifelse(acc/att > arate, 1.5, 0.7), 0.01), 10)
          tunep <- pmin(pmax(tunep*ifelse(acc/att > arate, 1.5, 0.7), 0.01), 10)
          att <- acc <- 0
      }

    ## SAVE

      if (iscan > nburn) {
          skipcount = skipcount + 1
          if (skipcount >= nskip) {
              isave <- isave + 1
              dispcount <- dispcount + 1
              gamma1save[isave, ] <- gamma1
              beta1save[isave] <- beta1
              sigma1save[isave] <- sigma1
              gamma2save[isave, ] <- gamma2
              betaysave[isave] <- betay
              sigma21save[isave] <- sigma21
              psave[isave] <- p
              skipcount <- 0
              if (dispcount >= ndisp) {
                  dispcount <- 0
                  cat(isave, 'out of', nsave, proc.time()[3] - start, '\n')
              }
          }
      }

  }

  ans <- list(gamma1save = gamma1save,
              beta1save = beta1save,
              sigma1save = sigma1save,
              gamma2save = gamma2save,
              betaysave = betaysave,
              sigma21save = sigma21save,
              psave = psave,
              n = n,
              xdim = xdim,
              y = y,
              X = X,
              R = R,
              tau = tau,
              tune = list(gamma1 = tunegamma1, beta1 = tunebeta1,
                sigma1 = tunesigma1, p = tunep, gamma2 =tunegamma2,
                betay = tunebetay, sigma21 = tunesigma21)
              )

  class(ans) <- 'QRMissingBiBayes'

  return(ans)

}

##' Observed Log Likelihood
##'
##'
##' @title observed log likelihood for bivariate case with single normal
##' @param gamma1
##' @param beta1
##' @param sigma1
##' @param gamma2
##' @param beta2sp
##' @param sigma21
##' @param sigma21sp
##' @param betay
##' @param betaysp
##' @param p
##' @param tau
##' @param y
##' @param X
##' @param R
##' @return log likelihood of bivariate response with missingness
##' @useDynLib qrmissing
##' @author Minzhao Liu, Mike Daniels
##' @export
LLBiSingle <- function(gamma1, beta1, sigma1,
                       gamma2, beta2sp, sigma21, sigma21sp, betay, betaysp,
                       p, tau, y, X, R){
  n <- dim(y)[1]
  xdim <- dim(X)[2]
  num <- sum(R)

  d <- matrix(0, n, 2)
  d <- .Fortran("mydelta2bise",
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

  lp1 <- beta1
  mu11 <- d[, 1] + lp1
  mu10 <- d[, 1] - lp1
  mu21 <- d[, 2] + betay * y[, 1]
  ll11 <- sum(dnorm(y[, 1], mu11, sigma1, log=T)[R==1])
  ll10 <- sum(dnorm(y[, 1], mu10, sigma1, log=T)[R==0])
  ll21 <- sum(dnorm(y[, 2], mu21, sigma21, log = T)[R==1])
  ans <- ll11 + ll10 + ll21 + num*log(p) + (n - num)*log(1 - p)

  return(ans)
}


##' @rdname QRMissingBiBayes
##' @method coef QRMissingBiBayes
##' @S3method coef QRMissingBiBayes
coef.QRMissingBiBayes <- function(mod, ...){
  gamma1 <- apply(mod$gamma1save, 2, mean)
  beta1 <- mean(mod$beta1save)
  sigma1 <- mean(mod$sigma1save)
  gamma2 <- apply(mod$gamma2save, 2, mean)
  betay <- mean(mod$betaysave)
  sigma21 <- mean(mod$sigma21save)
  p <- mean(mod$psave)
  return(list(gamma1 = gamma1, beta1 = beta1, sigma1 = sigma1,
              gamma2 = gamma2, betay = betay, sigma21 = sigma21, p = p))
}

##' @rdname QRMissingBiBayes
##' @method summary QRMissingBiBayes
##' @S3method summary QRMissingBiBayes
summary.QRMissingBiBayes <- function(mod, ...){
  n <- mod$n
  R <- mod$R
  tau <- mod$tau
  param <- mod$par
  q <- mod$xdim

  cat('Number of observations: ', n, '\n')
  cat('Sample proportion of observed data: ', sum(R)/n, '\n')
  cat('Estimated pi:', coef(mod)$p, '\n')
  cat('Quantile: ', tau, '\n')
  cat('Quantile regression coefficients: \n')
  print(coef(mod))
}

##' @rdname QRMissingBiBayes
##' @method plot QRMissingBiBayes
##' @S3method plot QRMissingBiBayes
plot.QRMissingBiBayes <- function(mod, ...){
  xdim <- mod$xdim
  for (i in 1:xdim){
    plot(ts(mod$gamma1save[, i]), main = paste('gamma1', i, sep = ''))
  }
  for (i in 1:xdim){
    plot(ts(mod$gamma2save[, i]), main = paste('gamma2', i, sep = ''))
  }
  for (i in 1){
    plot(ts(mod$beta1save), main = paste('beta1', i, sep = ''))
  }
  for (i in 1){
    plot(ts(mod$sigma1save), main = paste('sigma1', i, sep = ''))
  }
  plot(ts(mod$sigma21save), main = 'Sigma21')
  plot(ts(mod$psave), main = 'p')
}
