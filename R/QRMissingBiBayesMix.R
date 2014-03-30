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

    ## using M-H sampling method
    gamma1save <- gamma2save <- matrix(0, nsave, xdim)
    beta1save <- beta2spsave <- betaysave <- psave <- rep(0, nsave)
    mu1save <- mu2save <- sigma1save <- sigma2save <- omega1save <- omega2save <- matrix(0, nsave, K)

    ## TUNE
    tunegamma1 <- tunegamma2 <- rep(0.001, xdim)
    tunebeta2sp <- tunebetay <- tunebeta1 <- 0.001
    tunesigma1 <- tunesigma2 <- 0.001
    tunep <- 0.0003
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
    omega1 <- omega2 <- rep(1/K, K)
    betay <- betaysp <- 0
    p <- num/n
    G1 <- sample(K, size = n, replace = T)
    G2 <- sample(K, size = n, replace = T)

    isave <- 0
    skipcount <- 0
    dispcount <- 0
    nscan <- nburn + nskip * nsave

    start <- proc.time()[3]

########################################

    ## first
    loglikeo <- LLBiMix(gamma1, beta1, gamma2, beta2sp, mu1, sigma1,
                        mu2, sigma2, omega1, omega1, omega2, omega2,
                        betay, 0, p, tau, y, X, R, K, G1, G2)


    ## roll


    for (iscan in 1:nscan) {
        ## update new parameters
        loglikeo <- LLBiMix(gamma1, beta1, gamma2, beta2sp, mu1, sigma1, mu2, sigma2, omega1, omega1, omega2, omega2, betay, 0, p, tau, y, X, R, K, G1, G2)

        attgamma1 <- attgamma1 + 1
        gamma1c <- rnorm(xdim, gamma1, tunegamma1)
        beta1c <- rnorm(1, beta1, tunebeta1)
        gamma2c <- rnorm(xdim, gamma2, tunegamma2)
        beta2spc <- rnorm(1, beta2sp, tunebeta2sp)
        mu1c <- rnorm(K, mu1, 0.003)
        mu2c <- rnorm(K, mu2, 0.003)
        sigma1c <- pmax(0.01, rnorm(K, sigma1, 0.001))
        sigma2c <- pmax(0.01, rnorm(K, sigma2, 0.001))
        ## delta omega as small uniform shift
        dl <- 50
        da <- runif(1, min=-1,max=1)/dl
        db <- runif(1, min = -1, max=1)/dl
        domega1 <- c(da, db, -da - db)
        omega1c <- domega1 + omega1
        if (any(omega1c < 0) | any(omega1c > 1)) omega1c <- omega1
        da <- runif(1, min=-1,max=1)/dl
        db <- runif(1, min = -1, max=1)/dl
        domega2 <- c(da, db, -da - db)
        omega2c <- domega2 + omega2
        if (any(omega2c < 0) | any(omega2c > 1)) omega2c <- omega2
        betayc <- rnorm(1, betay, tunebetay)
        pc <- max(min(rnorm(1, p, tunep), 0.99), 0.01)

        ## ll of candidate
        loglikec <- LLBiMix(gamma1c, beta1c, gamma2c, beta2spc, mu1c, sigma1c, mu2c, sigma2c, omega1c, omega1c, omega2c, omega2c, betayc, 0, pc, tau, y, X, R, K, G1, G2)

        ## prior
        logpriorc <- sum(dnorm(gamma1c, gammapm, gammapv, log = T)) +
            dnorm(beta1c, betapm, betapv, log = T) +
                sum(dnorm(gamma2c, gammapm, gammapv, log = T)) +
                    dnorm(beta2spc, betapm, betapv, log = T) +
                        sum(dnorm(mu1c, mupm, mupv, log = T)) +
                            sum(dgamma(sigma1c, sigmaa, sigmab, log = T)) +
                                sum(dnorm(mu2c, mupm, mupv, log = T)) +
                                    sum(dgamma(sigma2c, sigmaa, sigmab, log = T)) +
                                        log(ddirichlet(omega1c, alpha)) +
                                            log(ddirichlet(omega2c, alpha)) +
                                                dnorm(betayc, betapm, betapv, log = T) +
                                                    dbeta(pc, alpha1, alpha2, log = T)


        logprioro <- sum(dnorm(gamma1, gammapm, gammapv, log = T)) +
            dnorm(beta1, betapm, betapv, log = T) +
                sum(dnorm(gamma2, gammapm, gammapv, log = T)) +
                    dnorm(beta2sp, betapm, betapv, log = T) +
                        sum(dnorm(mu1 , mupm, mupv, log = T)) +
                            sum(dgamma(sigma1 , sigmaa, sigmab, log = T)) +
                                sum(dnorm(mu2 , mupm, mupv, log = T)) +
                                    sum(dgamma(sigma2 , sigmaa, sigmab, log = T)) +
                                        log(ddirichlet(omega1 , alpha)) +
                                            log(ddirichlet(omega2 , alpha)) +
                                                dnorm(betay , betapm, betapv, log = T) +
                                                    dbeta(p , alpha1, alpha2, log = T)


        ## accept
        ratio <- loglikec + logpriorc - loglikeo - logprioro
        if (log(runif(1)) <= ratio) {
            accgamma1 <- accgamma1 + 1
            loglikeo <- loglikec
            gamma1 <- gamma1c
            beta1 <- beta1c
            gamma2 <- gamma2c
            beta2sp <- beta2spc
            mu1 <- mu1c
            sigma1 <- sigma1c
            mu2 <- mu2c
            sigma2 <- sigma2c
            omega1 <- omega1c
            omega2 <- omega2c
            betay <- betayc
            p <- pc
        }


        ## Update G1, G2,
        dd <- matrix(0, n, 2)
        dd <- .Fortran("mydelta2bisemix",
                       x = as.double(X),
                       gamma1 = as.double(gamma1),
                       beta1 = as.double(beta1),
                       sigma1 = as.double(sigma1),
                       gamma2 = as.double(gamma2),
                       beta2sp = as.double(beta2sp),
                       mu1 = as.double(mu1),
                       sigma1 = as.double(sigma1),
                       mu2 = as.double(mu2),
                       sigma2 = as.double(sigma2),
                       omega11 = as.double(omega1),
                       omega10 = as.double(omega1),
                       omega21 = as.double(omega2),
                       omega20sp = as.double(omega2),
                       betay = as.double(betay),
                       betaysp = as.double(0),
                       p = as.double(pi),
                       tau = as.double(tau),
                       n = as.integer(n),
                       xdim = as.integer(xdim),
                       delta = as.double(dd),
                       K = as.integer(K))$delta
        dd <- matrix(dd, n, 2)

        G1 <- G2 <- rep(0, n)
        for (i in 1:n) {
            if (R[i] == 1){
                prob1 <- omega1 * dnorm(y[i, 1], dd[i, 1] + beta1 + mu1, sigma1)
                prob2 <- omega2 * dnorm(y[i, 2], dd[i, 2] + beta2sp + betay * y[i, 1] + mu2, sigma2)
            } else {
                prob1 <- omega1 * dnorm(y[i, 1], dd[i, 1] - beta1 + mu1, sigma1)
                prob2 <- omega2
            }
            G1[i] <- rcat(1, prob1)
            G2[i] <- rcat(1, prob2)
        }


        ## TUNE

        ## SAVE

        if (iscan > nburn) {
            skipcount = skipcount + 1
            if (skipcount >= nskip) {
                isave <- isave + 1
                dispcount <- dispcount + 1
                gamma1save[isave, ] <- gamma1
                beta1save[isave] <- beta1
                gamma2save[isave, ] <- gamma2
                beta2spsave[isave] <- beta2sp
                betaysave[isave] <- betay
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
                beta2spsave = beta2spsave,
                sigma2save = sigma2save,
                psave = psave,
                n = n,
                xdim = xdim,
                y = y,
                X = X,
                R = R,
                K = K,
                tau = tau,
                tune = list(gamma1 = tunegamma1, beta1 = tunebeta1,
                    sigma1 = tunesigma1, p = tunep, gamma2 =tunegamma2,
                    betay = tunebetay)
                )

    class(ans) <- 'QRMissingBiBayesMix'

    return(ans)

}

##' @rdname QRMissingBiBayesMix
##' @method coef QRMissingBiBayesMix
##' @S3method coef QRMissingBiBayesMix
coef.QRMissingBiBayesMix <- function(mod, ...){
    gamma1 <- apply(mod$gamma1save, 2, mean)
    beta1 <- mean(mod$beta1save)
    gamma2 <- apply(mod$gamma2save, 2, mean)
    beta2sp <- mean(mod$beta2spsave)
    betay <- mean(mod$betaysave)
    p <- mean(mod$psave)
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
        plot(ts(mod$gamma1save[, i]), main = paste('gamma1', i, sep = ''))
    }
    for (i in 1:xdim){
        plot(ts(mod$gamma2save[, i]), main = paste('gamma2', i, sep = ''))
    }
    plot(ts(mod$psave), main = 'p')
}
