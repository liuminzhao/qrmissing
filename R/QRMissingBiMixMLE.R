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
QRMissingBiMixMLE <- function(formula, R, tau = 0.5, sp = NULL,
                              init = NULL, method = 'uobyqa',
                              control = list(maxit = 1000,
                                  trace = 0), hess = FALSE,
                              K = 2,
                              model = 'slope'){
    ## convert formula
    modfor <- model.frame(formula)
    y <- model.response(modfor, "numeric")
    X <- model.matrix(modfor)

    ## data
    n <- dim(y)[1]
    num <- sum(R)
    xdim <- dim(X)[2]

    ## initial
    if (is.null(sp)) {
        if (model == 'int'){
            sp <- 0
        } else if (model == 'slope'){
            sp <- rep(0, xdim)
        }
    }
    if (!is.null(init)){
        param <- init
    } else {
        lmcoef1 <- coef(rq(y[,1] ~ X[,-1], tau = tau))
        lmcoef2 <- coef(rq(y[,2][R == 1] ~ X[R == 1,-1], tau = tau))

        if (model == 'int' ) {
            param <- rep(0, 2*xdim + 3*K + 1)
            param[1:xdim] <- lmcoef1
            param[(xdim + 1):(2*xdim)] <- lmcoef2
            param[2*xdim + 2*K + 2] = qlogis(num/n)
        } else if (model == 'slope') {
            ## param = (q1(q), q2(q), sigma(K),
            ## omega(K-1), beta1(q), betay, p
            ## mu(K-1),
            ## thus dim(param) = 3q + 3K
            param <- rep(0, 3*xdim + 3*K)
            param[1:xdim] <- lmcoef1
            param[(xdim + 1):(2*xdim)] <- lmcoef2
            param[3*xdim + 3*K] = qlogis(num/n)
        }
    }

    ## nll
    nll <- function(param){
        ll2Mix(param, y, X, R, tau, sp, K, model)
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
    mod$model <- model

    class(mod) <- "QRMissingBiMixMLE"

    return(mod)

}

##' @rdname QRMissingBiMixMLE
##' @method coef QRMissingBiMixMLE
##' @S3method coef QRMissingBiMixMLE
coef.QRMissingBiMixMLE <- function(mod, ...){
    xdim <- q <- mod$xdim
    K <- mod$K
    model <- mod$model

    param <- mod$par
    gamma1 <- mod$par[c(1:q)]
    gamma2 <- mod$par[(q + 1):(2*q)]

    if (model == 'int') {
        sigma <- exp(param[(xdim*2 + 1):(xdim*2 + K)])
        omega <- rep(0, K)
        z <- plogis(param[(xdim*2 + K + 1):(xdim*2 + K*2 - 1)])
        z[K] <- 1
        omega <- z2omega(z)
        beta1 <- param[2*xdim + K*2]
        betay <- param[2*xdim + K*2 + 1] # for R = 1
        p <- plogis(param[2*xdim + K*2 + 2])
        mu <- param[(xdim*2 + 2*K + 3):(xdim*2 + K*3 + 1)]
    } else if (model == 'slope') {
        beta1 <- param[(2*q + 1):(3*q)]
        sigma <- exp(param[(3*q + 1):(q*3 + K)])
        omega <- rep(0, K)
        z <- plogis(param[(q*3 + K + 1):(q*3 + K*2 - 1)])
        z[K] <- 1
        omega <- z2omega(z)
        mu <- param[(xdim*3 + 2*K):(xdim*3 + K*3 - 2)]
        betay <- param[3*xdim + K*3 - 1]
        p <- plogis(param[3*xdim + K*3])
    }

    return(list(gamma1 = gamma1, gamma2 = gamma2,
                beta1 = beta1, betay = betay,
                sigma = sigma,
                omega = omega,
                p = p, mu = mu
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
    model <- mod$model
    if (model == 'int') {
        p <- plogis(param[2*xdim + K*2 + 2])
    } else if (model == 'slope') {
        p <- plogis(param[3*q + 3*K])
    }

    cat('Number of observations: ', n, '\n')
    cat('Sample proportion of observed data: ', sum(R)/n, '\n')
    cat('Estimated pi:', p, '\n')
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
ll2Mix <- function(param, y, X, R, tau, sp, K, model){
    n <- dim(y)[1]
    xdim <- q <- dim(X)[2]
    num <- sum(R)

    ## param = (q1(q), q2(q), sigma1(K), sigma2(K),
    ## omega1(K-1), omega2(K-1), beta1, betay, p),
    ## thus dim(param) = 2q + 4K + 1

    if (model == 'int'){
        gamma1 <- param[1:xdim]
        gamma2 <- param[(xdim + 1):(2*xdim)]
        sigma <- exp(param[(xdim*2 + 1):(xdim*2 + K)])
        omega <- rep(0, K)
        z <- plogis(param[(xdim*2 + K + 1):(xdim*2 + K*2 - 1)])
        z[K] <- 1
        omega <- z2omega(z)
        beta1 <- param[2*xdim + K*2]
        betay <- param[2*xdim + K*2 + 1] # for R = 1
        p <- plogis(param[2*xdim + K*2 + 2])
        mu <- param[(xdim*2 + 2*K + 3):(xdim*2 + K*3 + 1)]
    } else if (model == 'slope') {
        gamma1 <- param[1:xdim]
        gamma2 <- param[(xdim + 1):(2*xdim)]
        beta1 <- param[(2*q + 1):(3*q)]
        sigma <- exp(param[(3*q + 1):(q*3 + K)])
        omega <- rep(0, K)
        z <- plogis(param[(q*3 + K + 1):(q*3 + K*2 - 1)])
        z[K] <- 1
        omega <- z2omega(z)
        mu <- param[(xdim*3 + 2*K):(xdim*3 + K*3 - 2)]
        betay <- param[3*xdim + K*3 - 1]
        p <- plogis(param[3*xdim + K*3])
    }

    beta2sp <- sp # SP for R = 0
    sigma21sp <- 0  # SP for R = 0
    betaysp <- 0 # SP for R = 0

    mu[K] <- - sum(mu * omega[1:(K-1)])/omega[K]

    dd <- matrix(0, n, 2)

    if (model == 'int') {
        dd <- .Fortran("mydelta2bisemix",
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
                       betaysp = as.double(betaysp),
                       p = as.double(p),
                       tau = as.double(tau),
                       n = as.integer(n),
                       xdim = as.integer(xdim),
                       delta = as.double(dd),
                       K = as.integer(K))$delta
    } else if (model == 'slope') {
        dd <- .Fortran("mydelta2bisemixslope",
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
                       betaysp = as.double(betaysp),
                       p = as.double(p),
                       tau = as.double(tau),
                       n = as.integer(n),
                       xdim = as.integer(xdim),
                       delta = as.double(dd),
                       K = as.integer(K))$delta
    }

    dd <- matrix(dd, n, 2)

    if (model == 'int') {
        lp1 <- beta1
    } else if (model == 'slope') {
        lp1 <- X %*% beta1
    }

    mu11 <- dd[, 1] + lp1
    mu10 <- dd[, 1] - lp1
    mu21 <- dd[, 2] + betay * y[, 1]

    ans <- 0

    for (i in 1:n){
        if (R[i] == 1) {
            ans <- ans + log(p) + log(sum(omega * dnorm(y[i, 1], mu11[i] + mu, sigma))) + log(sum(omega * dnorm(y[i, 2], mu21[i] + mu, sigma)))
        } else {
            ans <- ans + log(1 - p) + log(sum(omega * dnorm(y[i, 1], mu10[i] + mu, sigma)))
        }
    }

    return(-ans)
}
