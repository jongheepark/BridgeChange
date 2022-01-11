## ---------------------------------------------------- #### ---------------------------------------------------- ################
## The idea is to develop a sparsity-induced prior-posterior model that fits data
## Bridge regression using mixture of normals representation.
## by JHP "Wed Oct 19 15:13:47 2016"
## ---------------------------------------------------- #### ---------------------------------------------------- ################
####################################
## here we goes! Main function....
####################################
#'
#' Hybrid Bridge Regression with Change-point
#'
#' Hybrid Univariate response linear change-point model with Bridge prior.
#'
#' @param y
#' Outcome vector.
#' @param X
#' Design matrix. Columns correspond to variables.
#' @param intercept
#' If \code{TRUE}, estimate intercept by MLE.
#' Default is \code{TRUE}.
#' This option does not affect the result when \code{n.break = 0}.
#' We recommend \code{intercept = TRUE} when the number of break is not zero.
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
#'
#' @param reduce.mcmc The number of reduced MCMC iterations for marginal likelihood computations.
#' If \code{reduce.mcmc = NULL}, \code{mcmc/thin} is used.

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
#' @param regime.duration.min The minimum length of time in each regime. Default is 5.
#' If regime is shorter than this limit, hybrid analysis is skipped.  
#' @param alpha.MH
#' If \code{TRUE}, alpha is updated by the Metropolis–Hastings algorithm.
#' If \code{FALSE} the Griddy gibbs sampler is used instead.
#' Default is \code{TRUE}.
#' @param ols.weight If TRUE, OLS estimates are used for adpative lasso weight vector. 
#' @param beta.alg
#' An algorithm to sample beta.
#' Default is \code{beta.alg = "SVD"}.
#' Also supported is \code{beta.alg = "BCK"} and \code{beta.alg = "CHL"}.
#' @param waic
#' If \code{TRUE}, WAIC is computed after the parameter estimation.
#' Default is \code{FALSE}.
#' @param marginal
#' If \code{TRUE}, the marginal likelihood is computed based on Chib's method.
#' Default is \code{FALSE}.
#' @return
#'
#' @name BridgeChangeRegHybrid
#' @importFrom mvtnorm rmvnorm
#' @importFrom copula retstable
#' @importFrom coda as.mcmc
#'
#' @useDynLib BridgeChange
#' @export

"BridgeChangeRegHybrid" <- function(y, X,                                             # inputs
                                    n.break = 0,                                      # number of breaks
                                    scale.data = TRUE,  intercept = TRUE,             # data transformations
                                    mcmc = 100, burn = 100, verbose = 100, thin = 1,  # mcmc related args
                                    reduce.mcmc = NULL, ols.weight = FALSE,
                                    c0 = 0.1, d0 = 0.1, nu.shape = 2.0, nu.rate = 2.0,# priors / hyper params
                                    known.alpha = FALSE, alpha.start = 1,             # alpha related args
                                    alpha.limit = FALSE, alpha.MH = TRUE,
                                    beta.start = NULL, beta.alg = "SVD",              # beta realted args
                                    regime.duration.min = 5,                          # minimum regime duration (lower limit)
                                    waic = FALSE, marginal = FALSE                    # model selection args
                                    ){
    
    ## ---------------------------------------------------- ##
    ##                preparing inputs                      ##
    ## ---------------------------------------------------- ##
    ## time stamp
    start.time <- proc.time(); call <- match.call()
    
    ## data transoformation
    X <- Xorig  <- as.matrix(X);
    resid <- NULL
    
    y.sd <- sd(y)
    X.sd <- apply(X, 2, sd)
    
    ## scaling first
    if (scale.data) {
        ydm <- scale(y)
        X   <- scale(X)
    }else{
        ydm <- y
    }
    
    ## common quantities
    K      <- ncol(X); ntime  <- length(y)
    m      <- n.break; ns     <- m + 1
    nstore <- mcmc / thin
    
    ## prior for transition matrix
    a <- NULL; b <- NULL
    
    ## ---------------------------------------------------- ##
    ##                   initialization                     ##
    ## ---------------------------------------------------- ##
    ## alpha0 <- runif(1)
    if(known.alpha){
        if(is.null(alpha.start)) {
            stop("When known.alpha = TRUE, a user must specify the value of alpha.start\n")
        } else {
            alpha = rep(alpha.start, ns)
        }
    }else{
        alpha.start <- 1
        alpha <- rep(alpha.start, ns)
    }
    
    lambda <- rmvnorm(ns, rep(1, K))
    tau    <- rep(1, ns)
    
    if (is.null(beta.start) & (K < ntime)) {
        ols    <-  lm(ydm ~ X - 1)
        ols.coef <- coef(ols);
        ols.vcov <- vcov(ols)
        if(sum(is.na(ols.coef)) > 0){
            cat("OLS regression of y on X generates NA's in coefficients! Random numbers are imputed for those NA's.")
            ols.coef[which(is.na(ols.coef))] <- rnorm(sum(is.na(ols.coef)))
            ols.vcov <- ifelse(is.na(ols.vcov), runif(1), ols.vcov)
        }
        
        ## set starting values
        beta   <-  rmvnorm(ns, ols.coef, ols.vcov)
        sig2   <-  rep(summary(ols)$sigma, ns)^2
    }else if (is.null(beta.start) & K >= ntime) {
        ## initialization for high dimensional case
        beta <- matrix(rnorm(K * ns), nrow = ns, ncol = K)
        sig2 <- 1 / rgamma(ns, 1, 1)
        
        ## development codes here -------------------------
                                        # cat("Initializing betas by SLOG\n")
                                        # if (intercept == TRUE) {
                                        #   beta_slog <- c(rnorm(1), SLOG(x = X, y = ydm, l = tau[1]))
                                        #
                                        # } else {
                                        #   beta_slog <- SLOG(x = X, y = ydm, l = tau[1])
                                        # }
                                        # beta <- matrix(beta_slog, nrow = ns, ncol = length(beta_slog), byrow = TRUE)
        
    }else{
        ## when starting values are assigned
        beta <- rmvnorm(ns, beta.start, diag(K))
        sig2 <- rep(1, ns)
    }
    
    if (n.break > 0) P <-  MCMCpack:::trans.mat.prior(m=m, n=ntime, a=0.9, b=0.1)   ## very quick change
    if (n.break > 0) A0 <- MCMCpack:::trans.mat.prior(m=m, n=ntime, a=a, b=b)
    
    ## ---------------------------------------------------- ##
    ##                    Setup for MCMC                    ##
    ## ---------------------------------------------------- ##
    ## holder
    alphadraws <- taudraws<- matrix(0, nstore, ns)
    sigmadraws <- beta0draws <- matrix(0, nstore, ns)
    betadraws  <- lambdadraws <- matrix(0, nstore, ns*K)
    Pmat       <- matrix(NA, nstore, ns)
    sdraws     <- matrix(data = 1, nstore, ntime)
    Z.loglike.array <- matrix(data = 0, nstore, ntime)
    
    ## prepare initial stuff
    ## ps.store <- matrix(0, T, ns)
    ps.store <- matrix(0, ntime, ns) ## ps.store <- rep(0, ntime)
    totiter <- mcmc + burn
    if (n.break > 0) {
        state <- sort(sample(1:ns, size = ntime, replace = TRUE, prob = (rep(1, ns))))
        ps <- matrix(NA, T, ns)
    } else {
        state <- rep(1, ntime)
        ps <- matrix(1, ntime, 1)
    }
    
    ## message with initial states
    if(verbose != 0) {
        cat("\n----------------------------------------------------\n")
        cat("MCMC Sampling of BridgeChangeReg Starts! \n")
        if (n.break > 0) cat("Initial state = ", table(state), "\n")
        ## cat("start.time: ", start.time, "\n")
        cat("----------------------------------------------------\n")
    }
    
    XX <- XY <- Xm <- Ym <- list()
    if(n.break == 0){
        nj  <-  ntime
        Xm[[1]] <- X; Ym[[1]] <- ydm
        XX[[1]] <- t(X) %*% X
        XY[[1]] <- t(X) %*% ydm
    }
    
    
    ## ---------------------------------------------------- ##
    ## Estimate intercept for initialization                ##
    ## ---------------------------------------------------- ##
    beta0 <- estimate_intercept_reg(y, Xorig, beta, n.break, intercept, state)
    
    ## ---------------------------------------------------- ##
    ## MCMC sampler starts here!
    ## ---------------------------------------------------- ##
    for (iter in 1:totiter) {
        ## if( i%%verbose==0 ) cat("iteration ", i, "\n")
        if(iter == (burn+1) ) {
            ess.time <- proc.time();
    }

    ## ---------------------------------------------------- ##
    ## Step 1: tau -- marginalized draw.
    ## ---------------------------------------------------- ##
    tau <- draw_tau_cpp(beta, alpha, nu.shape, nu.rate, ns)

    ## ---------------------------------------------------- ##
    ## Step 2: sig2
    ## ---------------------------------------------------- ##
    if (n.break > 0) {
      for (j in 1:ns){
        ej  <-  as.numeric(state==j)
        nj  <-  sum(ej)
        yj  <-  matrix(ydm[ej==1], nj, 1)
        Xj  <-  matrix(X[ej==1,], nj, K)
        sig2[j] <- draw.sig2(beta = beta[j,], x=Xj, y=yj, c0, d0)
        XX[[j]] <- crossprod(Xj, Xj)
        XY[[j]] <- crossprod(Xj, yj)
        Xm[[j]] <- Xj; Ym[[j]] <- yj
      }
    } else {
      sig2[1] <- draw.sig2(beta=beta[1,], x = X, y = ydm, c0, d0)
    }

    ## development code
    # sig2 <- draw_sig2_cpp(y, X, beta, state, c0, d0, ns)

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
          alpha[j] <- draw.alpha2(alpha[j], beta[j,], tau[j])
        }
      } else {
       ## if alpha is sampled from (0, 2]
        if (alpha.MH) {
          for (j in 1:ns){
            alpha[j] <- draw.alpha(alpha[j], beta[j,], tau[j]);
          }
        } else {
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
      if (intercept == TRUE){
        X1    <- cbind(rep(1, nrow(X)), X)
        beta1 <- cbind(beta0 / sd(y), beta) ## add back intercept to the design matrix
        state.out <- sparse_state_sampler_cpp(m, ydm, X1, beta1, sig2, P)
        state <- state.out$state
      } else {
        state.out <- sparse_state_sampler_cpp(m, ydm, X, beta, sig2, P)
        state <- state.out$state
      }
    }

    ## ---------------------------------------------------- ##
    ## Step 7: sampling P
    ## ---------------------------------------------------- ##
    if (n.break > 0) {
      switch <- switchg(state, m = m)
      ## cat("switch = ", print(switch) ,  "\n")
      for (j in 1:ns){
        switch1 <- A0[j,] + switch[j,]
        pj      <- rdirichlet.cp(1, switch1)
        P[j,]   <- pj
      }
    }


    ## ---------------------------------------------------- ##
    ## Estimate intercept
    ## ---------------------------------------------------- ##
    beta0 <- estimate_intercept_reg(y, Xorig, beta, n.break, intercept, state)

    ## ---------------------------------------------------- ##
    ## report and save
    ## ---------------------------------------------------- ##
    if (verbose != 0 & iter %% verbose == 0){
      cat(sprintf("\r Estimating parameters. Now at %i of %i", iter, totiter))
      flush.console()
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
      }

      if (waic) {
        mu.state <- X %*% t(beta)
        d <- sapply(1:ntime, function(t) {
          dnorm(ydm[t], mean = c(mu.state[t, state[t]]), sd = sqrt(sig2[state[t]]), log=TRUE)
        })
        Z.loglike.array[(iter-burn)/thin,] <- d
      }
    }

  }   ## end of MCMC iteration here

    ## ---------------------------------------------------- ##
    ## Marginal Likelihood Estimation starts here!
    ## ---------------------------------------------------- ##

    ## compute residual
    beta.st   <- matrix(apply(betadraws, 2, mean), ns, K, byrow=TRUE)
    mu.st.state <- X %*% t(beta.st)
    yhat.mat <- X%*%t(beta.st)
    prob.state <- cbind(sapply(1:ns, function(k){apply(sdraws == k, 2, mean)}))
    yhat <- apply(yhat.mat*prob.state, 1, sum)
    
    ## yhat <- sapply(1:ntime, function(t){c(mu.st.state[t, state[t]] + beta0[state[t]])})
    resid <- ydm - yhat
    Waic.out <- NULL
    if (isTRUE(marginal)) {

        ## run time information and waic
        end.time <- proc.time();
        runtime <- (end.time - start.time)[1];
        if (isTRUE(waic)){
            ## Waic computation
            Waic.out <- waic_calc(Z.loglike.array)$total
            rm(Z.loglike.array)

            if (verbose > 0) {
                cat("\n----------------------------------------------",'\n')
                cat("WAIC: ", Waic.out[1], "\n")
                ## cat("lpd: ", Waic.out[3], "\n")
                ## cat("p_Waic: ", Waic.out[4], "\n")
                cat("trun time: ", runtime, '\n')
                cat("----------------------------------------------",'\n')
            }
        } else {
            if (verbose > 0) {
                cat("\n----------------------------------------------",'\n')
                cat("trun time: ", runtime, '\n')
                cat("----------------------------------------------",'\n')
            }
            
        }
    }
           
    
    ## attr(output, "prob.state") <- ps.store/(mcmc/thin)
    ## pull together matrix and build MCMC object to return
    xnames <-  sapply(c(1:K), function(i) { paste("beta", i, sep = "") })
    lnames <- sapply(c(1:K), function(i) { paste("lambda", i, sep = "") })
    output <- NA
    
    if (n.break == 0) {
        if (isTRUE(scale.data)) betadraws <- betadraws * sd(y) / apply(X, 2, sd) ## report the coef in the original scale
        output <- as.mcmc(cbind(betadraws, sigmadraws))
        colnames(output) <- c(xnames, "sigma2")
    } else {
        sidx <- rep(1:ns, each = ncol(X))
        xidx <- 1:ncol(betadraws)
        idx  <- split(xidx, sidx)
            C1   <- y.sd / X.sd #sd(y) / apply(X, 2, sd)
        for (s in 1:ns) {
            betadraws[,idx[[s]]] <- t(apply(betadraws[,idx[[s]]], 1, function(x) x * C1))
        }
        
        output1 <- coda::mcmc(data=betadraws, start=burn+1, end=burn + mcmc, thin=thin)
        output2 <- coda::mcmc(data=sigmadraws, start=burn+1, end=burn + mcmc, thin=thin)
        lambda.out <- coda::mcmc(data=lambdadraws, start=burn+1, end=burn + mcmc, thin=thin)
        
        colnames(output1)  <- sapply(c(1:ns), function(i) { paste(xnames, "_regime", i, sep = "") })
        colnames(output2)  <- sapply(c(1:ns), function(i) { paste("sigma2_regime", i, sep = "") })
        colnames(lambda.out)  <- sapply(c(1:ns), function(i) { paste(lnames, "_regime", i, sep = "") })
        
        output    <- as.mcmc(cbind(output1, output2))
        ps.holder <- matrix(ps.store, ntime, ns)
        s.holder  <- matrix(sdraws, nstore, ntime)
    }

    ## ---------------------------------------------------- ##
    ## hybrid
    ## ---------------------------------------------------- ##
    ##  ninv.y <- 1/length(y)
    if(n.break > 0){
        state <- round(apply(sdraws, 2, mean))
        cat("\nEstiamted states are ", table(state), "\n")
        ## Following P. Richard HAHN and Carlos M. CARVALHO (Eq. 21)
        raw.y.list <- y.list <- x.list <- as.list(rep(NA, ns))
        ## raw.y.list.0 <- y.list.0 <- x.list.0 <- as.list(rep(NA, n.state))
        for(i in 1:ns){
            x.list[[i]] <- X[state==i, ]
            y.list[[i]] <- yhat[state==i]
            raw.y.list[[i]] <- y[state==i]   
            cat("y and yhat correlation is ", cor(y.list[[i]], raw.y.list[[i]]), "\n")
        }
        if(sum(lapply(y.list, length) < regime.duration.min )>0){

            if(ols.weight){
                hybrid.dss <- sapply(which(lapply(y.list, length)>=regime.duration.min),
                                    function(i){adaptive.lasso.olsweight(y.list[[i]], x.list[[i]])})
                rownames(hybrid.dss) <- colnames(X)
                colnames(hybrid.dss) <- paste0("Regime", which(lapply(y.list, length)>=regime.duration.min))
                
                hybrid.cp <- sapply(which(lapply(y.list, length)>=regime.duration.min),
                                  function(i){adaptive.lasso.olsweight(raw.y.list[[i]], x.list[[i]])})
                rownames(hybrid.cp) <- colnames(X)
                colnames(hybrid.cp) <- paste0("Regime", which(lapply(raw.y.list, length)>=regime.duration.min))
            }else{
                hybrid.dss <- sapply(which(lapply(y.list, length)>=regime.duration.min),
                                    function(i){adaptive.lasso(y.list[[i]], x.list[[i]], beta.hat = beta.st[i,])})
                rownames(hybrid.dss) <- colnames(X)
                colnames(hybrid.dss) <- paste0("Regime", which(lapply(y.list, length)>=regime.duration.min))
                
                hybrid.cp <- sapply(which(lapply(y.list, length)>=regime.duration.min),
                                  function(i){adaptive.lasso(raw.y.list[[i]], x.list[[i]], beta.hat = beta.st[i,])})
                rownames(hybrid.cp) <- colnames(X)
                colnames(hybrid.cp) <- paste0("Regime", which(lapply(raw.y.list, length)>=regime.duration.min))

            }
                
        }else{
            if(ols.weight){
                ## Variable selection using adaptive lasso
                hybrid.dss <- sapply(1:ns, function(i){adaptive.lasso.olsweight(y.list[[i]], x.list[[i]])})
                rownames(hybrid.dss) <- colnames(X)
                colnames(hybrid.dss) <- paste0("Regime", 1:ns)
                
                hybrid.cp <- sapply(1:ns, function(i){adaptive.lasso.olsweight(raw.y.list[[i]], x.list[[i]])})
                rownames(hybrid.cp) <- colnames(X)
                colnames(hybrid.cp) <- paste0("Regime", 1:ns)
            }else{
                ## Variable selection using adaptive lasso
                hybrid.dss <- sapply(1:ns, function(i){adaptive.lasso(y.list[[i]], x.list[[i]], beta.hat = beta.st[i,])})
                rownames(hybrid.dss) <- colnames(X)
                colnames(hybrid.dss) <- paste0("Regime", 1:ns)
                
                hybrid.cp <- sapply(1:ns, function(i){adaptive.lasso(raw.y.list[[i]], x.list[[i]], beta.hat = beta.st[i,])})
                rownames(hybrid.cp) <- colnames(X)
                colnames(hybrid.cp) <- paste0("Regime", 1:ns)
            }
        }
    }else{
        cat("y and yhat correlation is ", cor(yhat, y), "\n")
        if(ols.weight){
            hybrid.dss <- matrix(adaptive.lasso.olsweight(yhat, X), ncol(X), 1)
            rownames(hybrid.dss) <- colnames(X)
            
            hybrid.cp <- matrix(adaptive.lasso.olsweight(y, X), ncol(X), 1)
            rownames(hybrid.cp) <- colnames(X)
        }else{
            hybrid.dss <- matrix(adaptive.lasso(yhat, X, beta.hat = beta.st), ncol(X), 1)
            rownames(hybrid.dss) <- colnames(X)
            
            hybrid.cp <- matrix(adaptive.lasso(y, X, beta.hat = beta.st), ncol(X), 1)
            rownames(hybrid.cp) <- colnames(X)

        }
    }

    
    attr(output, "hybrid") <- hybrid.dss
    attr(output, "hybrid.raw") <- hybrid.cp 
    ## attr(output, "X") <- X
    attr(output, "title") <- "BridgeChangeReg Posterior Sample"
    attr(output, "intercept") <- coda::mcmc(beta0draws,start=burn+1, end=burn + mcmc, thin=thin)
    attr(output, "y")       <- y
    attr(output, "X")       <- X
    attr(output, "resid")   <- resid
    attr(output, "y.sd")    <- y.sd
    attr(output, "X.sd")    <- X.sd
    attr(output, "m")       <- m
    attr(output, "ntime")   <- ntime
    attr(output, "alpha")   <- coda::mcmc(data=alphadraws, start=burn+1, end=burn + mcmc, thin=thin)
    attr(output, "tau")     <- coda::mcmc(data=taudraws, start=burn+1, end=burn + mcmc, thin=thin)
    if (n.break > 0){
        attr(output, "s.store") <- sdraws
        prob.state <- cbind(sapply(1:ns, function(k){apply(sdraws == k, 2, mean)}))
        attr(output, "prob.state") <- prob.state
        attr(output, "lambda") <- lambda.out
    }
    attr(output, "Waic.out") <- Waic.out
    if(marginal){
        attr(output, "loglike") <- loglike
        attr(output, "logmarglike") <- logmarglike
    }
    
    class(output) <- c("mcmc", "BridgeChange")
    
    return(output)
} 
    
