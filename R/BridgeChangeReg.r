## ---------------------------------------------------- ##
## ---------------------------------------------------- ##
## The idea is to develop a sparsity-induced prior-posterior model that fits data
## Bridge regression using mixture of normals representation.
## by JHP "Wed Oct 19 15:13:47 2016"
## ---------------------------------------------------- ##
## ---------------------------------------------------- ##

#' Bridge Regression with Change-point
#'
#' Univariate response linear change-point model with Bridge prior.
#'
#' @param y
#' Outcome vector.
#' @param X
#' Design matrix. Columns correspond to variables.
#' @param intercept
#' If \code{TRUE}, estimate intercept by MLE.
#' Default is \code{TRUE}.
#' This option does not affect the result when \code{n.break = 0}.
#' We recommend \code{intercept = TRUE} when the number of break is not zero (change-point estimation).
#' @param n.break
#' The number of change-point(s).
#' If \code{n.break = 0}, the model corresponds to the usual regression.
#' Default is \code{n.break = 0}.
#' @param scale.data
#' If \code{TRUE}, \code{y} and \code{X} are both scaled to have zero mean and unit variance.
#' Default is \code{TRUE}.
#' We recommend \code{scale.data = TRUE} unless the original data are already scaled.
#' @param mcmc
#' The number of iterations for gibbs updates. Default is 100.
#' @param burn
#' The number of burn-in periods for gibbs updates. Default is 100.
#' @param verbose
#' Iterations at which the results are printed on the console. Default is every 100th iteration.
#' @param thin
#' Thinning for gibbs updates. Default is 1 (no thinning).
#' @param c0
#' Scale parameter for Gamma distribution. Used for the prior of \eqn{\sigma^2}.
#' Default is 0.1.
#' @param d0
#' Shape parameter for Gamma distribution. Used for the prior of \eqn{\sigma^2}.
#' Default is 0.1.
#' @param beta.start
#' Starting values of beta. If \code{NULL}, randomly choose beta.start from OLS or standard normal distribution.
#' Default is \code{NULL}.
#' @param nu.shape
#' Shape parameter for Gamma distribution. Used for the prior for \eqn{\tau}.
#' Default is 2.0.
#' @param nu.rate
#' Rate parameter for Gamma distribution. Used for the prior for \eqn{\tau}.
#' Default is 2.0.
#' @param alpha.limit
#' If \code{TRUE}, alpha is sampled from \eqn{(0,1)}, otherwise alpha is sampled between \eqn{(0,2)}.
#' Default is \code{FALSE}.
#' @param known.alpha
#' If \code{TRUE}, a user must specify a numeric value \eqn{[0, 2]} in \code{alpha.start}.
#'  Default is \code{FALSE} and therefore \eqn{\alpha} will be estimated.
#' @param alpha.start
#' Starting value for alpha.
#' When \code{known.alpha = TRUE}, alpha is fixed to the value of this argument.
#' Default is 1.
#' @param alpha.MH
#' If \code{TRUE}, alpha is updated by the Metropolis–Hastings algorithm.
#' If \code{FALSE} the Griddy gibbs sampler is used instead.
#' Default is \code{TRUE}.
#' @param beta.alg
#' An algorithm to sample beta.
#' Default is \code{beta.alg = "SVD"}.
#' Also supported is \code{beta.alg = "BCK"} and \code{beta.alg = "CHL"}.
#' @param Waic
#' If \code{TRUE}, WAIC is computed after the parameter estimation.
#' Default is \code{FALSE}.
#' @param marginal
#' If \code{TRUE}, the marginal likelihood is computed based on Chib's method.
#' Default is \code{FALSE}.
#' @return
#' An mcmc object that contains the posterior sample of coefficients and variances as columns.
#'  Rows correspond to mcmc samples. This object can be summarized by functions provided by the coda package.
#'  The object contains an attribute \code{"intercept"} that stores mcmc samples for intercept and an attribute state storage matrix that contains posterior samples of hidden states.
#' @example examples/reg_eg.R
#'
#' @importFrom mvtnorm rmvnorm
#' @importFrom copula retstable
#' @importFrom coda as.mcmc
#' @export
BridgeChangeReg <- function(
    y, X,                                             # inputs
    n.break = 0,                                      # break
    scale.data = TRUE,  intercept=TRUE,               # data transformations
    mcmc = 100, burn = 100, verbose = 100, thin = 1,  # mcmc related args
    c0 = 0.1, d0 = 0.1, nu.shape=2.0, nu.rate=2.0,    # priors
    beta.start = NULL,                                # initial values
    known.alpha = FALSE, alpha.start = 1,             # alpha related args
    alpha.limit = FALSE, alpha.MH = TRUE,
    beta.alg = "SVD",                                 # beta realted args
    Waic = FALSE, marginal = FALSE                    # model selection args
){
    start.time = proc.time();
    call <- match.call()
    ## require("copula"); ## No longer required.  Adapted implementation to BayesBridge.
    ## Set up.
    X <- Xorig  <- as.matrix(X);

    y.sd = sd(y)
    X.sd = apply(X, 2, sd)

    ## prior for transition matrix
    a = NULL; b = NULL

    # if(inter){
    #     x1.1 <- data.frame(X)
    #     var.names <- colnames(X)
    #     x1.2 <- matrix(t(apply(x1.1, 1, combn, 2, prod)), nrow = nrow(X))
    #     newX <- as.matrix(cbind(x1.1, x1.2))
    #     colnames(newX) <- c(var.names, combn(var.names, 2, paste, collapse="-"))
    #     X <- newX
    # }

    ## scaling first
    if (scale.data) {
        ydm <- scale(y)
        X   <- scale(X)
    }else{
        ydm <- y
    }

    # ## add constant for intercept estimation
    # if (n.break > 0 & intercept == TRUE) {
    #     X <- cbind(1, X)
    # }


    ## add constant for intercept estimation
    # if (intercept == TRUE) {
    #     X <- cbind(1, X)
    # }

    ## common quantities
    K      <- ncol(X)
    ntime  <- length(y)
    m      <- n.break
    ns     <- m + 1
    nstore <- mcmc/thin

    ## alpha0 <- runif(1)
    if(known.alpha){
        if(is.null(alpha.start)) {
            stop("When known.alpha = TRUE, a user must specify the value of alpha.start\n")
        } else {
            alpha = rep(alpha.start, ns)
        }
    }else{
        alpha.start = 1
        alpha = rep(alpha.start, ns);
    }

    ## Initialize
    lambda <- rmvnorm(ns, rep(1, K));
    tau    <- rep(1, ns);

    if(is.null(beta.start) & K < ntime){
        ols    <-  lm(ydm ~ X-1)
        beta   <-  rmvnorm(ns, coef(ols), vcov(ols))
        sig2   <-  rep(summary(ols)$sigma, ns)^2
    }else if(is.null(beta.start) & K >= ntime){
        beta <- matrix(rnorm(K * ns), nrow = ns, ncol = K)
        # cat("Initializing betas by SLOG\n")
        # if (intercept == TRUE) {
        #   beta_slog <- c(rnorm(1), SLOG(x = X, y = ydm, l = tau[1]))
        #
        # } else {
        #   beta_slog <- SLOG(x = X, y = ydm, l = tau[1])
        # }
        # beta <- matrix(beta_slog, nrow = ns, ncol = length(beta_slog), byrow = TRUE)
        sig2 <- 1 / rgamma(ns, 1, 1)
    }else{
        beta <- rmvnorm(ns, rnorm(K), diag(K))
        sig2 <- rep(1, ns)
    }

    if (n.break > 0) P <-  MCMCpack:::trans.mat.prior(m=m, n=ntime, a=0.9, b=0.1)   ## very quick change
    ## mvn.prior <- MCMCpack:::form.mvn.prior(b0, B0, K)
    ## b0 <- mvn.prior[[1]]
    ## B0 <- mvn.prior[[2]]

    ## prior
    if (n.break > 0) A0 <- MCMCpack:::trans.mat.prior(m=m, n=ntime, a=a, b=b)

    ## holder
    alphadraws  <- matrix(0, nstore, ns)
    taudraws    <- matrix(0, nstore, ns)
    betadraws   <- matrix(0, nstore, ns*K)
    lambdadraws <- matrix(0, nstore, ns*K)
    sigmadraws  <- matrix(0, nstore, ns)
    beta0draws  <- matrix(0, nstore, ns)
    Pmat        <- matrix(NA, nstore, ns)
    ## Ddraws <- matrix(data=0, nstore, ns*Q*Q)
    ## psdraws <- matrix(data=0, ntime, ns)
    sdraws <- matrix(data=0, nstore, ntime)
    Z.loglike.array <- matrix(data=0, nstore, ntime)

    ## make var names for each regime
    ## colnames(output$beta) = as.vector(sapply(1:ns, function(i){paste0(colnames(X), "_regime_", i)}))

    ## prepare initial stuff
    ## ps.store <- matrix(0, T, ns)
    ps.store <- matrix(0, ntime, ns) ## ps.store <- rep(0, ntime)
    totiter <- mcmc+burn
    if (n.break > 0) {
        state <- sort(sample(1:ns, size=ntime, replace=TRUE, prob=(rep(1, ns))))
        ## cat("\n randomly chosen initial state = ", table(state), "\n")
        ps <- matrix(NA, T, ns)
    } else {
        state <- rep(1, ntime)
        ps <- matrix(1, ntime, 1)
    }


    XX <- XY <- Xm <- Ym <- list()
    if(n.break==0){
        nj  <-  ntime
        Xm[[1]] <- X;
        Ym[[1]] <- ydm
        XX[[1]] <- t(X)%*%X
        XY[[1]] <- t(X)%*%ydm
    }

    if(verbose !=0){
        cat("\n----------------------------------------------------\n")
        cat("MCMC Sampling of BridgeChangeReg Starts! \n")
        cat("Initial state = ", table(state), "\n")
        ## cat("start.time: ", start.time, "\n")
        cat("----------------------------------------------------\n")
    }


    ## ---------------------------------------------------- ##
    ## Estimate intercept for initialization
    ## ---------------------------------------------------- ##
    if (n.break == 0 & intercept == TRUE) {
        ## fit intercept
        beta0 <- mean(y) - as.vector(colMeans(Xorig) %*% t(beta))
        # cat("beta0 = ", beta0, "\n")
    } else if (n.break > 0 & intercept == TRUE) {
        beta0 <- matrix(NA, nrow = ns, ncol = 1)
        # ydm   <- as.vector(as.vector(y) - tapply(y, state, mean)[state])
        for (j in 1:ns) {
            beta0[j,] <- mean(y[state == j]) - colMeans(Xorig[state == j, , drop= FALSE]) %*% beta[j,]
        }
    } else if (n.break > 0 & intercept == FALSE) {
        beta0 <- matrix(0, nrow = ns, ncol = 1)
        # ydm   <- as.vector(as.vector(y) - tapply(y, state, mean)[state])
    }

    ## ---------------------------------------------------- ##
    ## MCMC sampler starts here!
    ## ---------------------------------------------------- ##
    for(iter in 1:totiter){
        ## if( i%%verbose==0 ) cat("iteration ", i, "\n")
        if(iter == (burn+1) ) {
            ess.time = proc.time();
        }

        ## ---------------------------------------------------- ##
        ## Step 1: tau -- marginalized draw.
        ## ---------------------------------------------------- ##
        tau <- draw_tau_cpp(beta, alpha, nu.shape, nu.rate, ns)

        ## ---------------------------------------------------- ##
        ## Step 2: sig2
        ## ---------------------------------------------------- ##
        if(n.break > 0){
            for (j in 1:ns){
                ej  <-  as.numeric(state==j)
                nj  <-  sum(ej)
                yj  <-  matrix(ydm[ej==1], nj, 1)
                Xj  <-  matrix(X[ej==1,], nj, K)
                sig2[j] = draw.sig2(beta=beta[j,], x=Xj, y=yj, c0, d0)
                XX[[j]] <- crossprod(Xj, Xj)
                XY[[j]] <- crossprod(Xj, yj)
                Xm[[j]] <- Xj; Ym[[j]] <- yj
            }
        }else{
            sig2[1] = draw.sig2(beta=beta[1,], x=X, y=ydm, c0, d0)
        }
        ## for (j in 1:ns){
        # sig2[j] = draw.sig2(beta[j,], Xj, yj, c0, d0)
        # }
        ## cat("sig2[",j, "] = ", sig2[j], "\n");
        # sig2 <- draw_sig2_cpp(y, X, beta, state, c0, d0, ns)
        ## cat('sig2 = ', sig2, "\n")

        ## ---------------------------------------------------- ##
        ## Step 3: lambda
        ## ---------------------------------------------------- ##
        for (j in 1:ns){
            for(k in 1:K){
                ## lambda[j, k] =  2*retstable_LD(0.5 * alpha[j], 1.0, beta[j, k]^2 / tau[j]^2)
                ## lambda[j] = 2 * retstable.ld(0.5 * alpha, 1.0, beta[j]^2 / tau^2);
                lambda[j, k] =  2*retstable(0.5 * alpha[j], 1.0, beta[j, k]^2 / tau[j]^2, method="LD");
                # lambda[j, k] =  2*retstable_LD(0.5 * alpha[j], 1.0, beta[j, k]^2 / tau[j]^2)
                # cat("lambda = ", lambda[j,k], "\n")
            }
        }
        ## plot(1:K, lambda)

        ## ---------------------------------------------------- ##
        ## Step 4: beta
        ## ---------------------------------------------------- ##
        if(beta.alg %in% c("SVD")) {
            beta <- draw_beta_svd_cpp(Xm, Ym, lambda, sig2, tau, ns, K)
        } else if (beta.alg %in% c("BCK")) {
            beta <- draw_beta_BCK_cpp(Xm, Ym, lambda, sig2, tau, ns, K)
        } else if (beta.alg %in% c("CHL")){
            beta <- draw_beta_cpp(XX, XY, lambda, sig2, tau, ns, K)
        } else {
            stop("This algorithm is not supported. Please read the documentation on beta.alg!\n")
        }


        ## ---------------------------------------------------- ##
        ## Step 5: alpha
        ## ---------------------------------------------------- ##
        if (!known.alpha){
            if(alpha.limit){
                for (j in 1:ns){
                    # alpha[j] <- draw.alpha2(alpha[j], beta[j,], tau[j])
                    alpha[j] <- draw.alpha(alpha[j], beta[j,], tau[j], ub = 1.0)
                }
            }else{
                ## if alpha is sampled from (0, 2]
                if (alpha.MH) {
                    for (j in 1:ns){
                        alpha[j] <- draw.alpha(alpha[j], beta[j,], tau[j]);
                    }
                    ## alpha[j] <- draw.alpha.randomwalk(alpha[j], beta[j,], tau[j], window=0.1)
                }else {
                    for (j in 1:ns){
                        alpha[j] <- draw.alpha.griddy(beta[j,], tau[j])
                    }
                }
            }
        }

        ## ---------------------------------------------------- ##
        ## Step 6: sampling S
        ## ---------------------------------------------------- ##
        if (n.break > 0) {
            ## tau, alpha, lambda are all marginalized out in likelihood!
            ## state.out <- sparse.state.sampler(m=m, y=y, X=X, beta=beta, sig2=sig2, P = P)
            # if(sum(known.s)>0){
            #     state <- known.s
            #     ps    <- matrix(0, T, ns)
            # } else
            if (intercept == TRUE){
                X1    <- cbind(rep(1, nrow(X)), X)
                beta1 <- cbind(beta0 / sd(y), beta)
                state.out <- sparse_state_sampler_cpp(m, ydm, X1, beta1, sig2, P)
                state <- state.out$state
                ## if(length(table(state)) < ns){
                ##     state <- sort(sample(1:ns, size=ntime, replace=TRUE, prob=apply(ps, 2, mean)))
                ## }
                ## ps <- state.out$ps
            } else {
                state.out <- sparse_state_sampler_cpp(m, ydm, X, beta, sig2, P)
                state <- state.out$state
            }
        }
        ## table(state)

        ## ---------------------------------------------------- ##
        ## Step 7: sampling P
        ## ---------------------------------------------------- ##
        if (n.break > 0) {
            switch  <-  switchg(state, m=m)
            ## cat("switch = ", print(switch) ,  "\n")
            for (j in 1:ns){
                switch1 <-  A0[j,] + switch[j,]
                pj      <-  rdirichlet.cp(1, switch1)
                P[j,]   <-  pj
            }
        }

        ## ---------------------------------------------------- ##
        ## demeaning y by regime
        ## ---------------------------------------------------- ##
        # if (n.break > 0 & demean == TRUE) {
        #     ydm <- as.vector(as.vector(y) - tapply(y, state, mean)[state]) / tapply(y, state, sd)[state]
        # }

        ## ---------------------------------------------------- ##
        ## Estimate intercept
        ## ---------------------------------------------------- ##
        if (n.break == 0 & intercept == TRUE) {
            ## fit intercept
            beta0 <- mean(y)  - as.vector(colMeans(Xorig) %*% t(beta))
            # cat("beta  = ", beta, "\n")
            # cat("beta0 = ", beta0, "\n")
        } else if (n.break > 0 & intercept == TRUE) {
            beta0 <- matrix(NA, nrow = ns, ncol = 1)
            for (j in 1:ns) {
                beta0[j,] <- mean(y[state == j]) - colMeans(Xorig[state == j, , drop= FALSE]) %*% beta[j,]
            }

            # cat("beta0 = ", beta0, "\n")
        } else if (intercept == FALSE) {
            beta0 <- matrix(0, nrow = ns, ncol = 1)
            # ydm   <- as.vector(as.vector(y) - tapply(y, state, mean)[state])
        }

        ## ---------------------------------------------------- ##
        ## report
        ## ---------------------------------------------------- ##
       if (verbose!= 0 &iter %% verbose == 0){
            cat(sprintf("\r Estimating parameters. Now at %i of %i", iter, totiter))
            flush.console()
            # cat("\n----------------------------------------------",'\n')
            # cat("## iteration = ", iter, '\n')
            # cat("----------------------------------------------",'\n')
            # if (n.break > 0) {
            #     cat("sampled states: ", table(state), '\n')
            #     ## cat("Transition matrix: ", format(t(P), digits=2) , '\n')
            #     for(j in 1:ns){
            #         cat("alpha at state ", j, ": ", alpha[j], '\n')
            #         cat("beta at state ", j, ": ", format(beta[j, ], digits=2), '\n')
            #         cat("sigma^2 at state ", j, ": ", sig2[j] , '\n')
            #     }
            #
            # }else{
            #     cat("alpha: ", alpha, '\n')
            #     cat("beta: ", format(t(beta), digits=2), '\n')
            # }

        }
        if (iter > burn && (iter %% thin == 0)) {
            alphadraws[(iter-burn)/thin,]  <- alpha
            betadraws[(iter-burn)/thin,]   <- t(beta)
            lambdadraws[(iter-burn)/thin,] <- t(lambda)
            sigmadraws[(iter-burn)/thin,]  <- sig2
            taudraws[(iter-burn)/thin,]    <- tau
            beta0draws[(iter-burn)/thin,]  <- as.vector(beta0)
            if (n.break > 0) {
                Pmat[(iter-burn)/thin, ]   <- diag(P)
                sdraws[(iter-burn)/thin,]  <- state
                ## ps.store <- ps.store + ps
            }
            if(Waic){
                mu.state <- X %*% t(beta)
                ## py  <-  sapply(c(1:(m+1)), function(i){
                ## dnorm(y[t], mean= mu[i], sd = sqrt(sig2[i]))})

                d <- sapply(1:ntime, function(t){dnorm(ydm[t],
                                                      mean = c(mu.state[t, state[t]]),
                                                      sd=sqrt(sig2[state[t]]), log=TRUE)})
                Z.loglike.array[(iter-burn)/thin,] <- d
            }

        }

    }

    ## ---------------------------------------------------- ##
    ## Marginal Likelihood Estimation starts here!
    ## ---------------------------------------------------- ##
    if(marginal){

        ## ---------------------------------------------------- ##
        ## prepare
        ## ---------------------------------------------------- ##
        beta.st   <- matrix(apply(betadraws, 2, mean), ns, K, byrow=TRUE)
        lambda.st <- matrix(apply(lambdadraws, 2, mean), ns, K, byrow=TRUE)
        sig2.st   <- apply(sigmadraws, 2, mean)
        alpha.st  <- apply(alphadraws, 2, mean)
        tau.st    <- apply(taudraws, 2, mean)
        state.st  <- round(apply(sdraws, 2, median))

        if(n.break > 0){
            P.st <- apply(Pmat, 2, mean)
        }
        mu.st.state <- X %*% t(beta.st)

        ## ---------------------------------------------------- ##
        ## Likelihood computation
        ## ---------------------------------------------------- ##
        loglike.t <- sapply(1:ntime, function(t){dnorm(ydm[t],
                                                       mean = c(mu.st.state[t, state.st[t]]),
                                                       sd=sqrt(sig2.st[state.st[t]]), log=TRUE)})
        loglike <- sum(loglike.t)

        cat("\n---------------------------------------------- \n ")
        cat("Likelihood computation \n")
        cat("    loglike: ", as.numeric(loglike), "\n")
        cat("---------------------------------------------- \n ")

        ## ---------------------------------------------------- ##
        ## holders
        ## ---------------------------------------------------- ##
        density.sig2.holder <- density.nu.holder <- density.beta.holder <- density.lambda.holder <- matrix(NA, mcmc, ns)


        ## ---------------------------------------------------- ##
        ## Marginal Step 1. density.nu
        ## ---------------------------------------------------- ##
        ## draw.tau <- function(beta, alpha, c, d)
        ## {
        ##   p = length(beta)
        ##   nu = rgamma(1, c + p/alpha, rate=d + sum(abs(beta)^alpha))
        ##   tau = nu^(-1/alpha)
        ##   return(tau);
        ## }
        nu.st <- tau.st
        for(g in 1:mcmc){
            for(j in 1:ns){
                nu.st[j] <- tau.st[j]^(-alphadraws[g, j])
                shape.g  <- nu.shape + K/alphadraws[g, j]
                rate.g   <- nu.rate + sum(abs(betadraws[g,K*(j-1)+1:K])^alphadraws[g, j])
                density.nu.holder[g, j] <- dgamma(nu.st[j], shape.g, rate = rate.g)
            }
        }
        density.nu <- log(prod(apply(density.nu.holder, 2, mean)))
        cat("\n---------------------------------------------- \n ")
        cat("Marignal Likelihood Computation Step 1 \n")
        cat("    density.nu: ", as.numeric(density.nu), "\n")
        cat("---------------------------------------------- \n ")

        ## ---------------------------------------------------- ##
        ## Marginal Step 2: Sigma2|tau.st
        ## ---------------------------------------------------- ##
        for(g in 1:mcmc){
            ## Reduced Step 1: sig2
            for (j in 1:ns){
                ej  <-  as.numeric(state==j)
                nj  <-  sum(ej)
                yj  <-  matrix(ydm[ej==1], nj, 1)
                Xj  <-  matrix(X[ej==1,], nj, K)
                rss = sum( (as.matrix(yj)-Xj%*%as.matrix(beta[j,]))^2 )
                shape <- c0 + nj/2
                rate <- d0 + rss/2
                sig2[j] = rinvgamma(1, shape, scale=rate)
                density.sig2.holder[g, j] <- dinvgamma(sig2.st[j], shape, scale=rate)
                XX[[j]] <- t(Xj)%*%Xj
                XY[[j]] <- t(Xj)%*%yj
            }
            ## Fixed : tau
            ## tau <- draw_tau_cpp(beta, alpha, nu.shape, nu.rate, ns)

            ## Reduced Step 2: lambda (treating as latent)
            for (j in 1:ns){
                for(k in 1:K){
                    lambda[j, k] =  2*retstable(0.5 * alpha[j], 1.0, beta[j, k]^2 / tau.st[j]^2, method="LD");
                }
            }
            ## Reduced Step 3: beta
            beta <- draw_beta_svd_cpp(XX, XY, lambda, sig2, tau.st, ns, K)
            ## beta <- draw_beta_cpp(XX, XY, lambda, sig2, tau.st, ns, K)

            ## Reduced Step 4: alpha
            if (!known.alpha){
                if(alpha.limit){
                    for (j in 1:ns){
                        alpha[j] <- draw.alpha2(alpha[j], beta[j,], tau.st[j])
                    }
                }else{
                    ## if alpha is sampled from (0, 2]
                    if (alpha.MH) {
                        for (j in 1:ns){
                            alpha[j] <- draw.alpha(alpha[j], beta[j,], tau.st[j]);
                        }
                        ## alpha[j] <- draw.alpha.randomwalk(alpha[j], beta[j,], tau[j], window=0.1)
                    }else {
                        for (j in 1:ns){
                            alpha[j] <- draw.alpha.griddy(beta[j,], tau.st[j])
                        }
                    }
                }
            }

            ## Reduced Step 5: sampling S
            if(n.break > 0){
                state.out <- sparse_state_sampler_cpp(m, ydm, X, beta, sig2, P)
                state <- state.out$state
                ## if(length(table(state)) < ns){
                ##    state <- sort(sample(1:ns, size=ntime, replace=TRUE, prob=apply(ps, 2, mean)))
                ## }

            }

            ## Reduced Step 6: sampling P
            if(n.break > 0){
                switch  <-  switchg(state, m=m)
                for (j in 1:ns){
                    switch1 <-  A0[j,] + switch[j,]
                    pj      <-  rdirichlet.cp(1, switch1)
                    P[j,]  <-  pj
                }
            }

            # ## demeaning y by regime
            # if (n.break > 0 & demean == TRUE) {
            #     ydm <- as.vector(as.vector(y) - tapply(y, state, mean)[state])
            # }

        }
        density.sig2 <- log(prod(apply(density.sig2.holder, 2, mean)))
        cat("\n---------------------------------------------- \n ")
        cat("Marignal Likelihood Computation Step 2 \n")
        cat("    density.sig2: ", as.numeric(density.sig2), "\n")
        cat("---------------------------------------------- \n ")

        ## ---------------------------------------------------- ##
        ## Marginal Step 3: beta| tau.st, sig.st
        ## Note that we do not evaluate the posterior ordinate for beta
        ## as it is a latent variable with hyperparameter in Bayes Bridge model
        ## ---------------------------------------------------- ##
        ## for(g in 1:mcmc){
        ##      if(n.break > 0){
        ##         ## Fixed sig2
        ##         for (j in 1:ns){
        ##             ej  <-  as.numeric(state==j)
        ##             nj  <-  sum(ej)
        ##             yj  <-  matrix(ydm[ej==1], nj, 1)
        ##             Xj  <-  matrix(X[ej==1,], nj, K)
        ##             XX[[j]] <- t(Xj)%*%Xj
        ##             XY[[j]] <- t(Xj)%*%yj
        ##         }
        ##     }
        ## Fixed tau
        ## tau <- draw_tau_cpp(beta, alpha, nu.shape, nu.rate, ns)

        ## Reduced Step 1: lambda
        ##      for (j in 1:ns){
        ##          for(k in 1:K){
        ##              lambda[j, k] =  2*retstable(0.5 * alpha[j], 1.0, beta[j, k]^2 / tau.st[j]^2, method="LD");
        ##          }
        ##      }
        ## Reduced Step 2: beta
        ##      for (j in 1:ns){
        ##          VInv = (XX[[j]] + diag(lambda[j,] * sig2.st[j] / tau.st[j]^2, K));
        ##          V = chol2inv(chol(VInv));
        ##          U = chol(V) * sqrt(sig2.st[j]);
        ##          Mu = V %*% XY[[j]];
        ##          beta[j,] = drop(Mu + t(U) %*% rnorm(K))
        ##          density.beta.holder[g, j]  <- log(dmvnorm(beta.st[j,], Mu, V))
        ##     cat('beta [',j, ']', beta[j,], "\n")
        ##      }
        ## Reduced Step 3: alpha
        ##     for (j in 1:ns){
        ##         if (!alpha.MH) {
        ##              alpha[j] <- draw.alpha.griddy(beta[j,], tau.st[j])
        ##         } else {
        ##             alpha[j] <- draw.alpha(alpha[j], beta[j,], tau.st[j]);
        ##         }
        ##     }
        ## Reduced Step 4: sampling S
        ##     if(n.break > 0){
        ##         state.out <- sparse_state_sampler_cpp(m, ydm, X, beta, sig2.st, P)
        ##         state <- state.out$state
        ##         if(length(table(state)) < ns){
        ##             state <- sort(sample(1:ns, size=ntime, replace=TRUE, prob=apply(ps, 2, mean)))
        ##         }

        ##         ## Reduced Step  5: sampling P
        ##         switch  <-  switchg(state, m=m)
        ##         for (j in 1:ns){
        ##             switch1 <-  A0[j,] + switch[j,]
        ##             pj      <-  rdirichlet.cp(1, switch1)
        ##             P[j,]  <-  pj
        ##         }
        ##     }
        ## }
        ## density.beta <- sum(log(apply(exp(density.beta.holder), 2, mean)))
        ## cat("\n---------------------------------------------- \n ")
        ## cat("Marignal Likelihood Computation Step 3 \n")
        ## cat("    density.beta: ", as.numeric(density.beta), "\n")
        ## cat("---------------------------------------------- \n ")

        ## ---------------------------------------------------- ##
        ## Marginal Step 3: P| tau.st, sig.st
        ## ---------------------------------------------------- ##
        if(n.break > 0){
            density.P.holder <- matrix(NA, mcmc, ns-1)
            for(g in 1:mcmc){
                if(n.break > 0){
                    ## Fixed sig2
                    for (j in 1:ns){
                        ej  <-  as.numeric(state==j)
                        nj  <-  sum(ej)
                        yj  <-  matrix(ydm[ej==1], nj, 1)
                        Xj  <-  matrix(X[ej==1,], nj, K)
                        XX[[j]] <- t(Xj)%*%Xj
                        XY[[j]] <- t(Xj)%*%yj
                    }
                }

                ## Reduced Step 1: lambda
                for (j in 1:ns){
                    for(k in 1:K){
                        lambda[j, k] =  2*retstable(0.5 * alpha[j], 1.0, beta[j, k]^2 / tau.st[j]^2, method="LD");
                    }
                }
                ## Reduced Step 2: beta
                beta <- draw_beta_svd_cpp(XX, XY, lambda, sig2.st, tau.st, ns, K)

                ## beta <- draw_beta_cpp(XX, XY, lambda, sig2.st, tau.st, ns, K)

                ## Reduced Step 3: alpha
                if (!known.alpha){
                    if(alpha.limit){
                        for (j in 1:ns){
                            alpha[j] <- draw.alpha2(alpha[j], beta[j,], tau.st[j])
                        }
                    }else{
                        ## if alpha is sampled from (0, 2]
                        if (alpha.MH) {
                            for (j in 1:ns){
                                alpha[j] <- draw.alpha(alpha[j], beta[j,], tau.st[j]);
                            }
                            ## alpha[j] <- draw.alpha.randomwalk(alpha[j], beta[j,], tau[j], window=0.1)
                        }else {
                            for (j in 1:ns){
                                alpha[j] <- draw.alpha.griddy(beta[j,], tau.st[j])
                            }
                        }
                    }
                }
                ## Reduced Step 4: sampling S
                state.out <- sparse_state_sampler_cpp(m, y, X, beta, sig2.st, P)
                state <- state.out$state
                if(length(table(state)) < ns){
                    state <- sort(sample(1:ns, size=ntime, replace=TRUE, prob=apply(ps, 2, mean)))
                }


                ## Reduced Step  5: sampling P
                swit  <-  switchg(state, m=m)
                for (j in 1:ns){
                    swit1 <-  A0[j,] + swit[j,]
                    pj      <-  rdirichlet.cp(1, swit1)
                    P[j,]  <-  pj
                    if(j < ns){
                        shape1 <-  swit1[j]
                        shape2 <-  swit1[j + 1]
                        density.P.holder[g, j] <- dbeta(P.st[j], shape1, shape2)
                    }
                }

                ## demeaning y by regime
                ## if (n.break > 0 & demean == TRUE) {
                ##     ydm <- as.vector(as.vector(y) - tapply(y, state, mean)[state])
                ## }
            }

            density.P <- log(prod(apply(density.P.holder, 2, mean)))
            cat("\n---------------------------------------------- \n ")
            cat("Marignal Likelihood Computation Step 3 \n")
            cat("    density.P: ", as.numeric(density.P), "\n")
            cat("---------------------------------------------- \n ")
        }

        ## ---------------------------------------------------- ##
        ## prior density estimation
        ## ---------------------------------------------------- ##
        if(n.break > 0){
            expected.duration <- round(ntime/(m + 1))
            b <- 0.1
            a <- b * expected.duration
            density.P.prior <- rep(NA, ns-1)
            for (j in 1:ns){
                if(j < ns){
                    density.P.prior[j] <- log(dbeta(P.st[j], a, b)) ## p = 1
                }
            }
        }

        density.nu.prior <- density.sig2.prior <- rep(NA, ns)
        for (j in 1:ns){
            nu.st[[j]] <- tau.st[j]^(-alpha.st[j])
            density.nu.prior[j] <- log(dgamma(nu.st[[j]], nu.shape, rate = nu.rate))
            density.sig2.prior[j] <- log(dinvgamma(sig2.st[j], c0, d0))
        }

        ## ---------------------------------------------------- ##
        ## Log marginal likelihood addition
        ## ---------------------------------------------------- ##
        if(n.break > 0){
            logprior <- sum(density.nu.prior) + sum(density.sig2.prior) + sum(density.P.prior)
            logdenom <- (density.nu + density.sig2 + density.P)
        }else{
            logprior <- sum(density.nu.prior) + sum(density.sig2.prior)
            logdenom <- (density.nu + density.sig2)
        }
        logmarglike <- (loglike + logprior) - logdenom

        cat("\n----------------------------------------------------\n")
        cat("    log marginal likelihood = (loglike + logprior) - (density.parameters) \n")
        cat("    log marginal likelihood : ", as.numeric(logmarglike), "\n")
        cat("    log likelihood : ", as.numeric(loglike), "\n")
        cat("    logprior : ", as.numeric(logprior), "\n")
        cat("    log.nu.prior : ", as.numeric(sum(density.nu.prior)), "\n")
        cat("    log.sig2.prior : ", as.numeric(sum(density.sig2.prior)), "\n")
        if(n.break > 0){cat("    log.P.prior : ", as.numeric(sum(density.P.prior)), "\n")}
        cat("    log posterior density : ", as.numeric(logdenom), "\n")
        cat("----------------------------------------------------\n")
    }

    end.time = proc.time();
    ## end.time = proc.time();
    runtime = (end.time - start.time)[1];
    Waic.out <- NULL
    if(Waic == TRUE){
        ## Waic computation
        Waic.out <- waic(Z.loglike.array)$total
        rm(Z.loglike.array)

        cat("\n----------------------------------------------",'\n')
        cat("\tWaic: ", Waic.out[1], "\n")
        ## cat("lpd: ", Waic.out[3], "\n")
        ## cat("p_Waic: ", Waic.out[4], "\n")
        cat("\trun time: ", runtime, '\n')
        cat("----------------------------------------------",'\n')
    }else{
        cat("\n----------------------------------------------",'\n')
        cat("\trun time: ", runtime, '\n')
        cat("----------------------------------------------",'\n')

    }


    ## attr(output, "prob.state") <- ps.store/(mcmc/thin)
    ## pull together matrix and build MCMC object to return
    # if(intercept){
    #     xnames <-  sapply(c(1:(K-1)), function(i){paste("beta", i, sep = "")})
    #     lnames <- sapply(c(1:(K-1)), function(i){paste("lambda", i, sep = "")})
    # }else{
        xnames <-  sapply(c(1:K), function(i){paste("beta", i, sep = "")})
        lnames <- sapply(c(1:K), function(i){paste("lambda", i, sep = "")})
    # }
    output <- NA

    if (m == 0){
        if (scale.data) betadraws <- betadraws * sd(y) / apply(X, 2, sd)
        output <- as.mcmc(cbind(betadraws, sigmadraws))
        colnames(output) <- c(xnames, "sigma2")
    }else{
        sidx <- rep(1:ns, each = ncol(X))
        xidx <- 1:ncol(betadraws)
        idx  <- split(xidx, sidx)
        C1   <- y.sd / X.sd #sd(y) / apply(X, 2, sd)
        for (s in 1:ns) {
            betadraws[,idx[[s]]] <- t(apply(betadraws[,idx[[s]]], 1, function(x) x * C1))
        }

        output1 <- coda::mcmc(data=betadraws, start=burn+1, end=burn + mcmc, thin=thin)
        output2 <- coda::mcmc(data=sigmadraws, start=burn+1, end=burn + mcmc, thin=thin)
        output4 <- coda::mcmc(data=lambdadraws, start=burn+1, end=burn + mcmc, thin=thin)

        colnames(output1)  <- sapply(c(1:ns),
                                     function(i){
                                         paste(xnames, "_regime", i, sep = "")
                                     })
        colnames(output2)  <- sapply(c(1:ns),
                                     function(i){
                                         paste("sigma2_regime", i, sep = "")
                                     })
        colnames(output4)  <- sapply(c(1:ns),
                                     function(i){
                                         paste(lnames, "_regime", i, sep = "")
                                     })

        output    <- as.mcmc(cbind(output1, output2))
        ps.holder <- matrix(ps.store, ntime, ns)
        s.holder  <- matrix(sdraws, nstore, ntime)
    }
    
    ## attr(output, "X") <- X
    attr(output, "title") <- "BridgeChangeReg Posterior Sample"
    attr(output, "intercept") <- coda::mcmc(beta0draws,start=burn+1, end=burn + mcmc, thin=thin)
    attr(output, "y")       <- y
    attr(output, "X")       <- X
    attr(output, "y.sd")    <- y.sd
    attr(output, "X.sd")    <- X.sd
    attr(output, "m")       <- m
    attr(output, "ntime")   <- ntime
    attr(output, "alpha")   <- coda::mcmc(data=alphadraws, start=burn+1, end=burn + mcmc, thin=thin)
    attr(output, "tau")     <- coda::mcmc(data=taudraws, start=burn+1, end=burn + mcmc, thin=thin)
    if (m > 0){
        attr(output, "s.store") <- s.holder
        prob.state <- cbind(sapply(1:ns, function(k){apply(s.holder == k, 2, mean)}))
        attr(output, "prob.state") <- prob.state
        attr(output, "lambda") <- output4
    }
    attr(output, "Waic.out") <- Waic.out
    if(marginal){
        attr(output, "loglike") <- loglike
        attr(output, "logmarglike") <- logmarglike
    }

    class(output) <- c("mcmc", "BridgeChange")

    return(output)
}
