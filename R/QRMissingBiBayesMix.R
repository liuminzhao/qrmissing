##' Observed Log Likelihood
##'
##'
##' @title observed log likelihood for bivariate case with mixture normal
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
LLBiMix <- function(gamma1, beta1, gamma2, beta2sp, mu1, sigma1,
                    mu2, sigma2, omega11, omega10, omega21, omega20sp,
                    betay, betaysp, p, tau, y, X, R, K, G1, G2){
    n <- dim(y)[1]
    xdim <- dim(X)[2]
    num <- sum(R)

    dd <- matrix(0, n, 2)
    dd <- .Fortran("mydelta2bisemix",
                   x = as.double(X),
                   gamma1 = as.double(gamma1),
                   beta1 = as.double(beta1),
                   gamma2 = as.double(gamma2),
                   beta2sp = as.double(beta2sp),
                   mu1 = as.double(mu1),
                   sigma1 = as.double(sigma1),
                   mu2 = as.double(mu2),
                   sigma2 = as.double(sigma2),
                   omega11 = as.double(omega11),
                   omega10 = as.double(omega10),
                   omega21 = as.double(omega21),
                   omega20sp = as.double(omega20sp),
                   betay = as.double(betay),
                   betaysp = as.double(betaysp),
                   p = as.double(p),
                   tau = as.double(tau),
                   n = as.integer(n),
                   xdim = as.integer(xdim),
                   delta = as.double(dd),
                   K = as.integer(K))$delta
    dd <- matrix(dd, n, 2)

    lp1 <- beta1
    mu11 <- dd[, 1] + lp1
    mu10 <- dd[, 1] - lp1
    mu21 <- dd[, 2] + betay * y[, 1]
    ll11 <- sum(dnorm(y[, 1], mu11 + mu1[G1], sigma1[G1], log=T)[R==1])
    ll10 <- sum(dnorm(y[, 1], mu10 + mu1[G1], sigma1[G1], log=T)[R==0])
    ll21 <- sum(dnorm(y[, 2], mu21 + mu2[G2], sigma2[G2], log = T)[R==1])
    ans <- ll11 + ll10 + ll21 + num*log(p) + (n - num)*log(1 - p) +
        sum(log(omega11[G1])*(R == 1)) + sum(log(omega10[G1])*(R == 0)) +
            sum(log(omega21[G2])*(R==1)) + sum(log(omega20sp[G2])*(R==0))

    return(ans)
}

##' Function for quantile regression in the presence of monotone
##' missingness for Bivariate case
##'
##' Using Bayesian method for Bivariate quantile regression
##' with mixture normal error assumption
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
QRMissingBiBayesMix <- function(y, R, X, tau = 0.5,
                                mcmc, prior
                                ){

    ## data
    n <- dim(y)[1]
    xdim <- dim(X)[2]
    num <- sum(R)
    K <- 3 ## # of components

    ## prior
    gammapm <- prior$gammapm
    gammapv <- prior$gammapv
    betapm <- prior$betapm
    betapv <- prior$betapv
    mupm <- prior$mupm
    mupv <- prior$mupv
    sigmaa <- prior$sigmaa ## for sigma
    sigmab <- prior$sigmab
    alpha1 <- prior$alpha1 ## for p
    alpha2 <- prior$alpha2
    alpha <- prior$alpha ## for omega

    ## prior 1a

    ## prior 1b
    mupm <- 0
    mupv <- 1

    ## SP prior
    beta2pm <- prior$beta2pm
    beta2pv <- prior$beta2pv

    ## MCMC
    nsave <- mcmc$nsave
    nskip <- mcmc$nskip
    nburn <- mcmc$nburn
    ndisp <- mcmc$ndisp

    ## SAVE
    ## q = (gamma1(p), beta1, gamma2(p), beta2sp,  mu1(K),
    ## sigma1(K), mu2(K), sigma2(K), fomega1(K-1), fomega2(K-1),
    ## betay, fpi) ; dim(q) = 2p + 6K + 2
    ## external: G1, G2,
    ## Data: y, X, R, tau
    q <- rep(0, 2*xdim + K*6 + 2)
    qsave <- matrix(0, nsave, length(q))

    ## TUNE
    tunegamma1 <- tunebeta1 <- tunegamma2 <- rep(0.3, xdim)
    tunesigma1 <- tunesigma21 <- 0.3
    tunebetay <- 0.3
    tunep <- 0.1
    arate <- 0.25
    attgamma1 <- accgamma1 <- attgamma2 <- accgamma2 <- attbeta1 <- accbeta1 <- rep(0, xdim)
    attsigma1 <- attp <- accsigma1 <- accp <- attsigma21 <- accsigma21 <- 0
    attbetay <- accbetay <- 0

    ## initial
    gamma1 <- coef(rq(y[, 1] ~ X[, -1], tau = tau))
    beta1 <- beta2sp <- 0
    gamma2 <- coef(rq(y[R==1,2] ~ X[R==1, -1], tau = tau))
    mu1 <- rnorm(K)
    mu2 <- rnorm(K)
    sigma1 <- rep(1, K)
    sigma2 <- rep(1, K)
    omega1 <- omega2 <- rep(1, K-1)
    betay <- betaysp <- 0
    p <- num/n
    q <- c(gamma1, beta1, gamma2, beta2sp,
           mu1, sigma1, mu2, sigma2,
           omega1, omega2, betay, logit(p))
    G1 <- sample(K, size = n, replace = T)
    G2 <- sample(K, size = n, replace = T)

    isave <- 0
    skipcount <- 0
    dispcount <- 0
    nscan <- nburn + nskip * nsave

    ## initial for HMC
    epsilon <- 0.01
    L <- 15

    start <- proc.time()[3]

########################################

    ## remove 0
    y[R==0, 2] <- 0

    ## apply HMC
    qsave <- .Fortran("hmc",
                      y       = as.double(y),
                      X       = as.double(X),
                      R       = as.integer(R),
                      n       = as.integer(n),
                      xdim    = as.integer(xdim),
                      k       = as.integer(K),
                      tau     = as.double(tau),
                      nscan   = as.integer(nscan),
                      nburn   = as.integer(nburn),
                      ndisp   = as.integer(ndisp),
                      nskip   = as.integer(nskip),
                      nsave   = as.integer(nsave),
                      initial = as.double(q),
                      qsave   = as.double(qsave)
                      )$qsave

    qsave <- matrix(qsave, nrow = nsave, ncol = length(q))

    ## output
    ans <- list(qsave = qsave,
                n = n,
                xdim = xdim,
                y = y,
                X = X,
                R = R,
                K = K,
                tau = tau,
                tune = list(gamma1 = tunegamma1, beta1 = tunebeta1,
                    sigma1 = tunesigma1, p = tunep, gamma2 =tunegamma2,
                    betay = tunebetay, sigma21 = tunesigma21)
                )

    class(ans) <- 'QRMissingBiBayesMix'

    return(ans)

}

##' @rdname QRMissingBiBayesMix
##' @method coef QRMissingBiBayesMix
##' @S3method coef QRMissingBiBayesMix
coef.QRMissingBiBayesMix <- function(mod, ...){
    xdim <- mod$xdim
    K <- mod$K
    gamma1 <- apply(mod$qsave[, (1:xdim)], 2, mean)
    beta1 <- mean(mod$qsave[, (xdim + 1)])
    gamma2 <- apply(mod$qsave[, ((xdim + 2):(xdim*2 + 1))], 2, mean)
    beta2sp <- mean(mod$qsave[, (xdim*2 + 2)])
    betay <- mean(mod$qsave[, (xdim*2+K*6 + 1)])
    p <- inv.logit(mean(mod$qsave[, (xdim*2 + 2 + K*6)]))
    return(list(gamma1 = gamma1, beta1 = beta1, beta2sp = beta2sp,
                gamma2 = gamma2, betay = betay, p = p))
}

##' @rdname QRMissingBiBayesMix
##' @method summary QRMissingBiBayesMix
##' @S3method summary QRMissingBiBayesMix
summary.QRMissingBiBayesMix <- function(mod, ...){
  n <- mod$n
  R <- mod$R
  tau <- mod$tau
  param <- mod$par
  xdim <- mod$xdim

  cat('Number of observations: ', n, '\n')
  cat('Sample proportion of observed data: ', sum(R)/n, '\n')
  cat('Estimated pi:', coef(mod)$p, '\n')
  cat('Quantile: ', tau, '\n')
  cat('Quantile regression coefficients: \n')
  print(coef(mod))
}

##' @rdname QRMissingBiBayesMix
##' @method plot QRMissingBiBayesMix
##' @S3method plot QRMissingBiBayesMix
plot.QRMissingBiBayesMix <- function(mod, ...){
  xdim <- mod$xdim
  K <- mod$K
  for (i in 1:xdim){
    plot(ts(mod$qsave[, i]), main = paste('gamma1', i, sep = ''))
  }
  for (i in 1:xdim){
    plot(ts(mod$qsave[, (xdim*2 + 1 + i)]), main = paste('gamma2', i, sep = ''))
  }
  for (i in (xdim + 1)){
    plot(ts(mod$qsave[, i]), main = paste('beta1', i, sep = ''))
  }
  for (i in (xdim*2 + 2)){
    plot(ts(mod$qsave[, i]), main = paste('beta2sp', i, sep = ''))
  }
  for (i in (xdim*2+K*6 + 1)){
    plot(ts(mod$qsave[, i]), main = paste('betay', i, sep = ''))
  }

  plot(ts(inv.logit(mod$qsave[, xdim*2 + K*6 + 2])), main = 'p')
}
