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
##' @author Minzhao Liu, Mike Daniels
##' @export
QRMissingBiMixMLE <- function(y, R, X, tau = 0.5, sp = NULL,
                              init = NULL, method = 'uobyqa',
                              control = list(maxit = 1000,
                                  trace = 0), hess = FALSE,
                              K = 1){
    ## data
    n <- dim(y)[1]
    num <- sum(R)
    xdim <- dim(X)[2]

    ## initial
    if (is.null(sp)) {
        sp <- 0 # beta2sp
    }
    if (!is.null(init)){
        param <- init
    } else {
        lmcoef1 <- coef(rq(y[,1] ~ X[,-1], tau = tau))
        lmcoef2 <- coef(rq(y[,2][R == 1] ~ X[R == 1,-1], tau = tau))

        ## param = (q1(q), q2(q), sigma1(K), sigma2(K),
        ## omega1(K-1), omega2(K-1), beta1, betay, p
        ## mu1(K-1), mu2(K-1)),
        ## thus dim(param) = 2q + 4K + 1
        param <- rep(0, 2*xdim + 6*K - 1)
        param[1:xdim] <- lmcoef1
        param[(xdim + 1):(2*xdim)] <- lmcoef2
        param[2*xdim + 4*K + 1] = qlogis(num/n)
    }

    ## nll
    nll <- function(param){
        ll2Mix(param, y, X, R, tau, sp, K)
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
    ## res <- residuals(mod$par, y, X, R, tau, sp)
    res <- NULL

    ## Hessian matrix and grad
    if (hess) {
        Hessian <- hessian(nll, mod$par)
        d <- grad(nll, mod$par)
        Jninv <- solve(Hessian)
        se <- matrix(0, 2, xdim)
        se[1, ] <- sqrt(diag(Jninv)[1:xdim])
        se[2, ] <- sqrt(diag(Jninv)[(xdim + 1):(2*xdim)])
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
    mod$K <- K

    class(mod) <- "QRMissingBiMixMLE"

    return(mod)

}

##' @rdname QRMissingBiMixMLE
##' @method coef QRMissingBiMixMLE
##' @S3method coef QRMissingBiMixMLE
coef.QRMissingBiMixMLE <- function(mod, ...){
    xdim <- q <- mod$xdim
    K <- mod$K
    param <- mod$par
    gamma1 <- mod$par[c(1:q)]
    gamma2 <- mod$par[(q + 1):(2*q)]
    sigma1 <- exp(param[(xdim*2 + 1):(xdim*2 + K)])
    sigma2 <- exp(param[(xdim*2 + 1 + K):(xdim*2 + K*2)])
    omega1 <- omega2 <- rep(0, K)
    omega1[1:(K-1)] <- exp(param[(xdim*2 + K*2 + 1):(xdim*2 + K*3 - 1)])/(sum(exp(param[(xdim*2 + K*2 + 1):(xdim*2 + K*3 - 1)])) + 1)
    omega1[K] <- 1/(sum(exp(param[(xdim*2 + K*2 + 1):(xdim*2 + K*3 - 1)])) + 1)
    omega2[1:(K-1)] <- exp(param[(xdim*2 + K*3):(xdim*2 + K*4 - 2)])/(sum(exp(param[(xdim*2 + K*3):(xdim*2 + K*4 - 2)])) + 1)
    omega2[K] <- 1/(sum(exp(param[(xdim*2 + K*3):(xdim*2 + K*4 - 2)])) + 1)
    beta1 <- param[2*xdim + K*4 - 1]
    betay <- param[2*xdim + K*4] # for R = 1
    p <- plogis(param[2*xdim + K*4 + 1])
    mu1 <- param[(xdim*2 + 4*K + 2):(xdim*2 + K*5)]
    mu2 <- param[(xdim*2 + 5*K + 1):(xdim*2 + K*6 - 1)]

    return(list(gamma1 = gamma1, gamma2 = gamma2,
                beta1 = beta1, betay = betay,
                sigma1 = sigma1, sigma2 = sigma2,
                omega1 = omega1, omega2 = omega2,
                p = p, mu1 = mu1, mu2 = mu2
                ))

}

##' @rdname QRMissingBiMixMLE
##' @method print QRMissingBiMixMLE
##' @S3method print QRMissingBiMixMLE
print.QRMissingBiMixMLE <- function(mod, ...){
    cat('Coefficients: \n')
    print(coef(mod))
}

##' @rdname QRMissingBiMixMLE
##' @method summary QRMissingBiMixMLE
##' @S3method summary QRMissingBiMixMLE
summary.QRMissingBiMixMLE <- function(mod, ...){
    n <- mod$n
    R <- mod$R
    tau <- mod$tau
    param <- mod$par
    q <- mod$xdim
    K <- mod$K

    cat('Number of observations: ', n, '\n')
    cat('Sample proportion of observed data: ', sum(R)/n, '\n')
    cat('Estimated pi:', plogis(param[2*q + 4*K + 1]), '\n')
    cat('Quantile: ', tau, '\n')
    cat('Number of components: ', K, '\n')
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
##' @rdname QRMissingBiMixMLE
##' @method plot QRMissingBiMixMLE
##' @S3method plot QRMissingBiMixMLE
plot.QRMissingBiMixMLE <- function(mod, ...){
    R <- mod$R
    res <- mod$res
    tau <- mod$tau
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
##' @export
ll2Mix <- function(param, y, X, R, tau, sp, K){
    n <- dim(y)[1]
    xdim <- dim(X)[2]
    num <- sum(R)

    ## param = (q1(q), q2(q), sigma1(K), sigma2(K),
    ## omega1(K-1), omega2(K-1), beta1, betay, p),
    ## thus dim(param) = 2q + 4K + 1
    gamma1 <- param[1:xdim]
    gamma2 <- param[(xdim + 1):(2*xdim)]
    sigma1 <- exp(param[(xdim*2 + 1):(xdim*2 + K)])
    sigma2 <- exp(param[(xdim*2 + 1 + K):(xdim*2 + K*2)])
    omega1 <- omega2 <- rep(0, K)
    omega1[1:(K-1)] <- exp(param[(xdim*2 + K*2 + 1):(xdim*2 + K*3 - 1)])/(sum(exp(param[(xdim*2 + K*2 + 1):(xdim*2 + K*3 - 1)])) + 1)
    omega1[K] <- 1/(sum(exp(param[(xdim*2 + K*2 + 1):(xdim*2 + K*3 - 1)])) + 1)
    omega2[1:(K-1)] <- exp(param[(xdim*2 + K*3):(xdim*2 + K*4 - 2)])/(sum(exp(param[(xdim*2 + K*3):(xdim*2 + K*4 - 2)])) + 1)
    omega2[K] <- 1/(sum(exp(param[(xdim*2 + K*3):(xdim*2 + K*4 - 2)])) + 1)
    beta1 <- param[2*xdim + K*4 - 1]
    betay <- param[2*xdim + K*4] # for R = 1
    p <- plogis(param[2*xdim + K*4 + 1])
    mu1 <- param[(xdim*2 + 4*K + 2):(xdim*2 + K*5)]
    mu2 <- param[(xdim*2 + 5*K + 1):(xdim*2 + K*6 - 1)]

    beta2sp <- sp # SP for R = 0
    sigma21sp <- 0  # SP for R = 0
    betaysp <- 0 # SP for R = 0

    mu1[K] <- - sum(mu1 * omega1[1:(K-1)])/omega1[K]
    mu2[K] <- - sum(mu2 * omega2[1:(K-1)])/omega2[K]

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
                   omega11 = as.double(omega1),
                   omega10 = as.double(omega1),
                   omega21 = as.double(omega2),
                   omega20sp = as.double(omega2),
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

    ans <- 0

    for (i in 1:n){
        if (R[i] == 1) {
            ans <- ans + log(p) + log(sum(omega1 * dnorm(y[i, 1], mu11[i] + mu1, sigma1))) + log(sum(omega2 * dnorm(y[i, 2], mu21[i] + mu2, sigma2)))
        } else {
            ans <- ans + log(1 - p) + log(sum(omega1 * dnorm(y[i, 1], mu10[i] + mu1, sigma1)))
        }
    }

    return(-ans)
}
