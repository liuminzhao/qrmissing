##' Function for quantile regression in the presence of monotone
##' missingness for univariate case
##'
##' Using Bayesian method for univariate case
##' for mixture of normals
##'
##' @title Quantile Regression in the Presence of Monotone Missingness; Univariate case with mixture normals; using HMC
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
QRMissingUniBayesMix <- function(y, R, X, tau = 0.5,
                                 mcmc, prior
                                 ){

    ## data
    n <- length(y)
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

    ## MCMC
    nsave <- mcmc$nsave
    nskip <- mcmc$nskip
    nburn <- mcmc$nburn
    ndisp <- mcmc$ndisp

    ## SAVE and initials
    ## betasave <- gammasave <- matrix(0, nsave, xdim)
    ## psave <- rep(0, nsave)
    ## sigmasave <- matrix(0, nsave, 2)
    q <- rep(0, xdim + K*3 + 2)
    qsave <- matrix(0, length(q), nsave)

    ## TUNE
    tunegamma <- tunebeta <- 0.3
    tunep <- 0.1
    tunesigma <- 0.3
    arate <- 0.25
    att <- acc <- 0

    ## initial

    gamma <- rep(0, xdim)
    beta <- 0
    mu <- rnorm(K)
    sigma <- rep(1, K)
    omega <- rep(1, K)
    p <- num/n
    q <- c(gamma, beta, mu, sigma, omega, logit(p))
    G <- sample(K, size = n, replace = T)

    isave <- 0
    skipcount <- 0
    dispcount <- 0
    nscan <- nburn + nskip * nsave

    ## initial for HMC
    epsilon <- 0.01
    L <- 15

    start <- proc.time()[3]

########################################

    ## roll

    ## q = (gamma, beta, mu, sigma, omega, p), G

    l <- length(q)

    U <- function(q){
        l <- length(q)
        xdim <- dim(X)[2]
        n <- length(y)

        gamma <- q[1:xdim]
        beta <- q[xdim + 1]
        mu <- q[(xdim + 2):(xdim + K + 1)]
        sigma <- exp(q[(xdim + K + 2):(xdim + K*2 + 1)])
        omega <- exp(q[(xdim + K*2 + 2):(xdim + K*3 +1)])/sum(exp(q[(xdim + K*2 + 2):(xdim + K*3 + 1)]))
        ## omega0 <- exp(q[(2*xdim + K*3 + 1):(2*xdim + K*4)])/sum(exp(q[(2*xdim + K*3 + 1):(2*xdim + K*4)]))
        pi <- exp(q[xdim + K*3 + 2])/(1 + exp(q[xdim + K*3 + 2]))

        ll <- LLUniMix(gamma, beta, K, mu, sigma, omega, omega,
                       pi, tau, y, X, R, G)

        prior1 <- sum(dnorm(gamma, gammapm, gammapv, log = T))
        prior2 <- sum(dnorm(beta, betapm, betapv, log = T))
        prior3 <- sum(dnorm(mu, log = T))
        prior4 <- 0 # sigma
        prior5 <- log(ddirichlet(omega, alpha))
        prior6 <- dbeta(pi, alpha1/2, alpha2/2, log = T)

        ll <- ll + prior1 + prior2 + prior3 + prior4 + prior5 + prior6

        return(-ll)

    }

    dU <- function(q){
        e <- 0.01
        oldU <- U(q)
        l <- length(q)
        oldU <- U(q)
        ans <- rep(0, l)

        for (i in 1:l){
            qprime <- q
            qprime[i] <- qprime[i] + e
            ans[i] <- (U(qprime) - oldU)/e
        }

        return(ans)

    }


    for (iscan in 1:nscan) {

        ## start HMC
        ## update q and p
        newq <- q
        p <- rnorm(l, 0, 1)
        current_p <- p

        p <- p - epsilon * dU(q)/2

        for (i in 1:L){
            newq <- newq + epsilon * p
            if (i!=L) p <- p - epsilon * dU(newq)
        }

        p <- p - epsilon * dU(newq)/2

        p <- -p

        current_U <- U(q)
        current_K <- sum(current_p^2)/2
        proposed_U <- U(newq)
        proposed_K <- sum(p^2/2)

        if (runif(1) < exp(current_U - proposed_U + current_K - proposed_K)){
            q <- newq
        }

        ## update G
        gamma <- q[1:xdim]
        beta <- q[xdim + 1]
        mu <- q[(xdim + 2):(xdim + K + 1)]
        sigma <- exp(q[(xdim + K + 2):(xdim + K*2 + 1)])
        omega <- exp(q[(xdim + K*2 + 2):(xdim + K*3 +1)])/sum(exp(q[(xdim + K*2 + 2):(xdim + K*3 + 1)]))
        ## omega0 <- exp(q[(2*xdim + K*3 + 1):(2*xdim + K*4)])/sum(exp(q[(2*xdim + K*3 + 1):(2*xdim + K*4)]))
        pi <- exp(q[xdim + K*3 + 2])/(1 + exp(q[xdim + K*3 + 2]))

        dd <- rep(0, n)
        dd <- .Fortran("mydelta1bisemix",
                       x = as.double(X),
                       gamma = as.double(gamma),
                       beta = as.double(beta),
                       K = as.integer(K),
                       mu = as.double(mu),
                       sigma = as.double(sigma),
                       omega1 = as.double(omega),
                       omega0 = as.double(omega),
                       pi = as.double(pi),
                       tau = as.double(tau),
                       n = as.integer(n),
                       xdim = as.integer(xdim),
                       delta = as.double(dd))$delta

        G <- rep(0, n)
        for (i in 1:n){
            if (R[i] == 1){
                prob <- omega * dnorm(y[i], dd[i] + beta + mu, sigma)
            } else {
                prob <- omega * dnorm(y[i], dd[i] - beta + mu, sigma)
            }
            G[i] <- rcat(1, prob)
        }

        ## SAVE

        if (iscan > nburn) {
            skipcount <- skipcount + 1
            if (skipcount >= nskip) {
                isave <- isave + 1
                dispcount <- dispcount + 1
                qsave[, isave] <- q
                skipcount <- 0
                if (dispcount >= ndisp) {
                    dispcount <- 0
                    cat(isave, 'out of', nsave, proc.time()[3] - start, '\n')
                }

            }
        }
    }

    ans <- list(qsave = qsave,
                n = n,
                xdim = xdim,
                y = y,
                X = X,
                R = R,
                tau = tau,
                K = K,
                tune = list(gamma = tunegamma, beta = tunebeta,
                    sigma = tunesigma, p = tunep)
                )

    class(ans) <- 'QRMissingUniBayesMix'

    return(ans)

}

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title G
##' @param ndraws
##' @param prob
##' @return G
##' @author Minzhao Liu
rcat <- function(ndraws, prob){
    ## generate categorical variables:
    prob[is.na(prob)]<-0
    if(sum(prob)==0){prob[1]<-1}
    (1:length(prob)) %*% rmultinom(ndraws, 1, prob)
}


##' Observed Log Likelihood
##'
##'
##' @title observed log likelihood for univariate case with mixture of normals
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
LLUniMix <- function(gamma, beta, K, mu, sigma,
                     omega1, omega0, p, tau, y, X, R, G){
    n <- length(y)
    xdim <- dim(X)[2]
    num <- sum(R)

    d <- rep(0, n)
    d <- .Fortran("mydelta1bisemix",
                  x = as.double(X),
                  gamma = as.double(gamma),
                  beta = as.double(beta),
                  K = as.integer(K),
                  mu = as.double(mu),
                  sigma = as.double(sigma),
                  omega1 = as.double(omega1),
                  omega0 = as.double(omega0),
                  p = as.double(p),
                  tau = as.double(tau),
                  n = as.integer(n),
                  xdim = as.integer(xdim),
                  delta = as.double(d))$delta

    lp <- beta
    mu11 <- d + lp
    mu10 <- d - lp
    ll11 <- sum(dnorm(y, mu11 + mu[G], sigma[G], log=T)[R==1])
    ll10 <- sum(dnorm(y, mu10 + mu[G], sigma[G], log=T)[R==0])
    ans <- ll11 + ll10 + num*log(p) + (n - num)*log(1 - p) + sum(log(omega1[G])*(R == 1)) +
        sum(log(omega0[G])*(R == 0))

    return(ans)
}




##' @rdname QRMissingUniBayesMix
##' @method coef QRMissingUniBayesMix
##' @S3method coef QRMissingUniBayesMix
coef.QRMissingUniBayesMix <- function(mod, ...){
    xdim <- mod$xdim
    K <- mod$K
    gamma <- apply(mod$qsave[(1:xdim), ], 1, mean)
    beta <- mean(mod$qsave[(xdim + 1), ])
    p <- inv.logit(mean(mod$qsave[(xdim + 2 + K*3), ]))
    return(list(gamma = gamma, beta = beta, p = p))
}

##' @rdname QRMissingUniBayesMix
##' @method summary QRMissingUniBayesMix
##' @S3method summary QRMissingUniBayesMix
summary.QRMissingUniBayesMix <- function(mod, ...){
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

##' @rdname QRMissingUniBayesMix
##' @method plot QRMissingUniBayesMix
##' @S3method plot QRMissingUniBayesMix
plot.QRMissingUniBayesMix <- function(mod, ...){
    xdim <- mod$xdim
    K <- mod$K
    for (i in 1:xdim){
        plot(ts(mod$qsave[i, ]), main = paste('gamma', i, sep = ''))
    }
    for (i in (xdim + 1)){
        plot(ts(mod$qsave[i, ]), main = paste('beta', i, sep = ''))
    }
    plot(ts(inv.logit(mod$qsave[xdim + 2 + K*3, ])), main = 'p')
}
