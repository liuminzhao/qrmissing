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
                                mcmc, prior, method = "DP"
                                ){

    ## data
    n <- dim(y)[1]
    xdim <- dim(X)[2]
    num <- sum(R)

    ## replace NA with 0
    y[is.na(y)] <- 0

    ## prior for \xi
    gammapm <- prior$gammapm
    gammapv <- prior$gammapv
    betapm <- prior$betapm
    betapv <- prior$betapv
    alpha1 <- prior$alpha1 ## for p
    alpha2 <- prior$alpha2
    ## TODO alpha prior
    ## alpha <- prior$alpha ## for z
    K <- prior$K


    ## constant
    sigmamu <- prior$sigmamu
    sigmaa <- prior$sigmaa ## for sigma ## nu1; gamma scale
    sigmab <- prior$sigmab ## nu2 ; rate
    eta1 <- prior$eta1 ## alpha prior
    eta2 <- prior$eta2
    A <- prior$A

    ## SP prior
    beta2pm <- prior$beta2pm
    beta2pv <- prior$beta2pv

    ## MCMC
    nsave <- mcmc$nsave
    nskip <- mcmc$nskip
    nburn <- mcmc$nburn
    ndisp <- mcmc$ndisp

    ## SAVE
    ## q = (gamma1(p), beta1, gamma2(p), beta2sp,  mu(K-1),
    ## sigma(K), z(K-1),
    ## betay, fpi) ; dim(q) = 2p + 3K + 2
    ## external: G1, G2,
    ## Data: y, X, R, tau

    ## using M-H sampling method
    gamma1save <- gamma2save <- matrix(0, nsave, xdim)
    beta1save <- beta2spsave <- betaysave <- psave <- rep(0, nsave)
    musave <- sigmasave <- omegasave <- matrix(0, nsave, K)
    alphasave <- thetasave <- rep(0, nsave)

    ## TUNE
    tunegamma1 <- prior$tunegamma1
    tunegamma2 <- prior$tunegamma2
    tunebeta2sp <- prior$tunebeta2sp
    tunebetay <- prior$tunebetay
    tunebeta1 <- prior$tunebeta1
    tunep <- prior$tunep
    arate <- prior$arate
    tunemu <- prior$tunemu
    tuneomega <- prior$tuneomega
    tunesigma <- prior$tunesigma

    attgamma1 <- accgamma1 <- attgamma2 <- accgamma2 <- attbeta1 <- accbeta1 <- rep(0, xdim)
    attsigma <- attp <- accsigma <- accp <- attsigma21 <- accsigma21 <- 0
    attbetay <- accbetay <- 0
    att <- acc <- 0

    ## initial
    gamma1 <- coef(rq(y[, 1] ~ X[, -1], tau = tau))
    beta1 <- beta2sp <- 0
    gamma2 <- coef(rq(y[R==1,2] ~ X[R==1, -1], tau = tau))
    mu <- rep(0, K)
    ## mu <- sort(mu, decreasing = TRUE)
    sigma <- rep(1, K)
    omega <- rep(1/K, K)
    z <- 1/(K:1)
    betay <- betaysp <- 0
    p <- num/n
    G1 <- sample(K, size = n, replace = T)
    G2 <- sample(K, size = n, replace = T)

    alpha <- 1
    theta <- 0

    ## constraint
    tol <- 0.01
    mu[K] <- - sum((mu * omega)[-K])/omega[K]
    if (omega[K] < tol) mu[K] <- 0

    dd <- matrix(0, n, 2)

    isave <- 0
    skipcount <- 0
    dispcount <- 0
    nscan <- nburn + nskip * nsave

    start <- proc.time()[3]

########################################

    ## first
    loglikeo <- LLBiMix(gamma1, beta1, gamma2, beta2sp, mu, sigma, mu, sigma, omega, omega, omega, omega, betay, 0, p, tau, y, X, R, K, G1, G2)

    ## roll

    for (iscan in 1:nscan) {
        ## update new parameters
        loglikeo <- LLBiMix(gamma1, beta1, gamma2, beta2sp, mu, sigma, mu, sigma, omega, omega, omega, omega, betay, 0, p, tau, y, X, R, K, G1, G2)

        att <- att + 1
        gamma1c <- rnorm(xdim, gamma1, tunegamma1)
        beta1c <- rnorm(1, beta1, tunebeta1)
        gamma2c <- rnorm(xdim, gamma2, tunegamma2)
        beta2spc <- beta2sp
        muc <- rnorm(K, mu, tunemu)
        if (method == "scale") {
            muc <- rep(0, K)
        }
        ## muc <- sort(muc, decreasing = TRUE)
        sigmac <- pmax(0.01, rnorm(K, sigma, tunesigma))
        ## delta omega as small uniform shift
        dl <- tuneomega
        dz <- runif(K, min = -1, max=1)/dl
        dz[K] <- 0
        zc <- z + dz
        if (any(zc < 0) | any(zc > 1)) zc <- z
        omegac <- z2omega(zc)
        ## TODO mu constraint
        muc[K] <- - sum((muc * omegac)[-K])/omegac[K]
        if (omegac[K] < tol) muc[K] <- 0
        if (method == "scale") {
            muc <- rep(0, K)
        }
        betayc <- rnorm(1, betay, tunebetay)
        pc <- max(min(rnorm(1, p, tunep), 0.99), 0.01)

        ## ll of candidate
        loglikec <- LLBiMix(gamma1c, beta1c, gamma2c, beta2spc, muc, sigmac, muc, sigmac, omegac, omegac, omegac, omegac, betayc, 0, pc, tau, y, X, R, K, G1, G2)

        ## prior
        logpriorc <- sum(dnorm(gamma1c, gammapm, gammapv, log = T)) +
            dnorm(beta1c, betapm, betapv, log = T) +
                sum(dnorm(gamma2c, gammapm, gammapv, log = T)) +
                    dnorm(beta2spc, betapm, betapv, log = T) +
                        sum(dnorm(muc[-K], theta, sigmamu, log = T)) +
                            sum(dgamma(sigmac^-1, sigmaa, rate = sigmab, log = T)) +
                                sum(dbeta(zc[-K], 1, alpha, log = T)) +
                                    dnorm(betayc, betapm, betapv, log = T) +
                                        dbeta(pc, alpha1, alpha2, log = T)

        logprioro <- sum(dnorm(gamma1, gammapm, gammapv, log = T)) +
            dnorm(beta1, betapm, betapv, log = T) +
                sum(dnorm(gamma2, gammapm, gammapv, log = T)) +
                    dnorm(beta2sp, betapm, betapv, log = T) +
                        sum(dnorm(mu[-K], theta, sigmamu, log = T)) +
                            sum(dgamma(sigma^-1 , sigmaa, scale = sigmab, log = T)) +
                                sum(dbeta(z[-K], 1, alpha, log = T)) +
                                    dnorm(betay , betapm, betapv, log = T) +
                                        dbeta(p , alpha1, alpha2, log = T)

        ## accept
        ratio <- loglikec + logpriorc - loglikeo - logprioro
        if (log(runif(1)) <= ratio) {
            acc <- acc + 1
            loglikeo <- loglikec

            gamma1 <- gamma1c
            beta1 <- beta1c
            gamma2 <- gamma2c
            beta2sp <- beta2spc
            mu <- muc
            sigma <- sigmac
            z <- zc
            omega <- omegac
            betay <- betayc
            p <- pc
        }

        ## update beta2sp separately
        beta2sp <- rnorm(1, beta2pm, beta2pv)

        ## Update G1, G2,
        updateg1g2 <- .Fortran("updateg",
                               x = as.double(X),
                               gamma1 = as.double(gamma1),
                               beta1 = as.double(beta1),
                               gamma2 = as.double(gamma2),
                               beta2sp = as.double(beta2sp),
                               mu1 = as.double(mu),
                               sigma1 = as.double(sigma),
                               mu2 = as.double(mu),
                               sigma2 = as.double(sigma),
                               omega11 = as.double(omega),
                               omega10 = as.double(omega),
                               omega21 = as.double(omega),
                               omega20sp = as.double(omega),
                               betay = as.double(betay),
                               betaysp = as.double(0),
                               p = as.double(pi),
                               tau = as.double(tau),
                               n = as.integer(n),
                               xdim = as.integer(xdim),
                               delta = as.double(dd),
                               K = as.integer(K),
                               y = as.double(y),
                               r = as.integer(R),
                               g1 = as.integer(G1),
                               g2 = as.integer(G2))


        G1 <- updateg1g2$g1
        G2 <- updateg1g2$g2

        ## update alpha
        alpha <- rgamma(1, K + eta1 - 1, rate = eta2 - sum(log(1 - z[-K])))

        ## update theta
        sigmastar <- sqrt(1/((K-1)/sigmamu^2 + 1/A^2))
        thetastar <- sigmastar^2 * sum(mu[-K])/sigmamu^2
        theta <- rnorm(1, thetastar, sigmastar)

        ## TUNE
        ## if (att >= 100 && iscan < nburn) {
        ## No tune
        if (FALSE) {
            prop <- acc/att
            tunegamma1 <- tunegamma1*ifelse(prop > arate, 1.5, 0.75)
            tunebeta1 <- tunebeta1*ifelse(prop > arate, 1.5, 0.75)
            tunegamma2 <- tunegamma2*ifelse(prop > arate, 1.5, 0.75)
            tunebeta2sp <- tunebeta2sp*ifelse(prop > arate, 1.5, 0.75)
            tunebetay <- tunebetay*ifelse(prop > arate, 1.5, 0.75)
            tunep <- tunep*ifelse(prop > arate, 1.5, 0.75)
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
                gamma2save[isave, ] <- gamma2
                beta2spsave[isave] <- beta2sp
                betaysave[isave] <- betay
                musave[isave,] <- mu
                sigmasave[isave, ] <- sigma
                omegasave[isave, ] <- omega
                psave[isave] <- p
                alphasave[isave] <- alpha
                thetasave[isave] <- theta
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
                gamma2save = gamma2save,
                betaysave = betaysave,
                beta2spsave = beta2spsave,
                musave = musave,
                sigmasave = sigmasave,
                omegasave = omegasave,
                alphasave = alphasave,
                thetasave = thetasave,
                psave = psave,
                n = n,
                xdim = xdim,
                y = y,
                X = X,
                R = R,
                K = K,
                tau = tau,
                tune = list(gamma1 = tunegamma1, beta1 = tunebeta1, sigma1 = tunesigma, p = tunep, gamma2 =tunegamma2, betay = tunebetay),
                mcmc = mcmc
                )

    class(ans) <- 'QRMissingBiBayesMix'

    return(ans)

}

##' @rdname QRMissingBiBayesMix
##' @method coef QRMissingBiBayesMix
##' @S3method coef QRMissingBiBayesMix
coef.QRMissingBiBayesMix <- function(mod, ...){
    nsave <- mod$mcmc$nsave
    nburn <- mod$mcmc$nburn

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
plot.QRMissingBiBayesMix <- function(mod, full = FALSE, ...){
    xdim <- mod$xdim
    K <- mod$K
    for (i in 1:xdim){
        plot(ts(mod$gamma1save[, i]), main = paste('gamma1', i, sep = ''))
    }
    for (i in 1:xdim){
        plot(ts(mod$gamma2save[, i]), main = paste('gamma2', i, sep = ''))
    }

    plot(ts(mod$beta2spsave), main = 'beat2sp')
    plot(ts(mod$betaysave), main = 'beaty')
    plot(ts(mod$psave), main = 'p')

    if (full) {
        for (i in 1:K){
            plot(ts(mod$musave[, i]), main = 'mu')
        }

        for (i in 1:K){
            plot(ts(mod$sigmasave[, i]), main = 'sigma')
        }

        for (i in 1:K){
            plot(ts(mod$omegasave[, i]), main = 'omega')
        }
    }

}
##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title convert z to omega
##' @param z
##' @return omega
##' @author Minzhao Liu
z2omega <- function(z){
    K <- length(z)
    omega <- rep(0, K)
    omega[1] <- z[1]
    for (i in 2:K){
        omega[i] <- z[i]*prod(1 - z[1:(i - 1)])
    }
    return(omega)
}
