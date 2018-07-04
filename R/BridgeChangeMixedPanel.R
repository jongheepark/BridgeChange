## ---------------------------------------------------- #### ---------------------------------------------------- ################
## The idea is to develop a sparsity-induced prior-posterior model that fits data
## Bridge regression using mixture of normals representation.
## by JHP "Wed Oct 19 15:13:47 2016"
## ---------------------------------------------------- #### ---------------------------------------------------- ################
####################################
## here we goes! Main function....
####################################

#' Mixed Effect Panel Model with Change-point
#'
#' Linear mixed effect model with change-points.
#'
#' @param subject.id unit id.
#' @param time.id time id.
#' @param y Reponse vector
#' @param X Design matrix for fixed effects
#' @param W Design matrix for random effects.
#' @param n.break Number of structural breaks (changepoints). If \code{n.break = 0} is specified,
#'  usual panel model is run with shrinkage prior on coefficients.
#' @param mcmc =100
#' @param burn =100
#' @param verbose =100 Verbose
#' @param thin Thinning
#' @param b0 Hyperparam
#' @param B0 Hyperparam
#' @param c0 = 0.1
#' @param d0 = 0.1
#' @param r0 Hyperparam
#' @param R0 Hyperparam
#' @param a = NULL
#' @param b = NULL
#' @param nu.shape =2.0
#' @param nu.rate =2.0
#' @param alpha = 1
#' @param alpha.MH If \code{TRUE}, alpha is updated by MH algorithm. By default, it is updated by griddy gibbs.
#' @param standardize Default is \code{TRUE}.
#'  If \code{TRUE}, Y will be demeaned, and X and W will be scaled to be mean zero and variance one.
#'  Estimated parameters are always in the transformed scale.
#' @param beta.start \code{= NA}
#' @param sigma2.start = NA
#' @param D.start= NA
#' @param P.start = NA
#' @param Waic \code{= TRUE} If \code{= TRUE}, compute WAIC and store as "Waic.out." Retrieve using attr().
#' @param marginal \code{= FALSE}
#' @param fixed \code{= TRUE} If \code{= TRUE}, the fixed effects model is chosen as a baseline model.
#'
#' @return output
#'
#' @importFrom MCMCpack rinvgamma rwish
#' @importFrom mvtnorm rmvnorm
#' @importFrom copula retstable
#' @importFrom coda as.mcmc mcmc
#' @useDynLib BridgeChange
#' @export
BridgeMixedPanel <- function(
    y, X, W,
    subject.id,
    time.id,
    standardize = TRUE,
    n.break = 1,
    mcmc=100, burn=100, verbose=100, thin = 1,
    b0, B0, c0 = 0.1, d0 = 0.1, r0, R0, a = NULL, b = NULL,
    nu.shape=2.0, nu.rate=2.0, alpha = 1, alpha.MH = FALSE,
    beta.start = NULL, sigma2.start = NA, D.start= NA, P.start = NA,
    Waic = FALSE, marginal = FALSE, fixed = TRUE,
    unscaled.Y = unscaled.Y, unscaled.X = unscaled.X
) {
    ## ---------------------------------------------------- ##
    # Data
    ## ---------------------------------------------------- ##
    call <- match.call()
    m  <- n.break;
    ns <- m + 1
    NT <-  length(y)
    K  <-  ncol(X)
    Q  <-  ncol(W)
    N  <- length(unique(subject.id))  # number of subject
    T  <- length(unique(time.id))## length(y)

    ## ---------------------------------------------------- ##
    # data transformation
    ## ---------------------------------------------------- ##
    if (standardize) {
        ## save original information 
        ysd <- sd(y)
        Xsd <- apply(X, 2, sd)
        
        ## demeaning Y
        Y <- scale(y) #as.matrix(y - mean(y, na.rm = TRUE))
        X <- as.matrix(scale(X))
        if (!fixed) W <- as.matrix(scale(W))
    } else {
        Y <- as.matrix(y)
        X <- as.matrix(X)
        W <- as.matrix(W)
    }

    ## ---------------------------------------------------- ##
    # take all pair-wise interactions
    ## ---------------------------------------------------- ##
    # if(inter){
    #     if(length(which(apply(X, 2, sd)==0)) > 0){
    #         X0 <- X[, -which(apply(X, 2, sd)==0)]
    #     }else{
    #         X0 <- X
    #     }
    #     x1.1 <- data.frame(X0)
    #     var.names <- colnames(X0)
    #     x1.2 <- t(apply(x1.1, 1, combn, 2, prod))
    #     newX <- as.matrix(cbind(x1.1, x1.2))
    #     colnames(newX) <- c(var.names, combn(var.names, 2, paste, collapse="-"))
    #     X <- cbind(X[, which(apply(X, 2, sd)==0)], newX)
    #     K  <-  ncol(X)
    # 
    # }


    ## ---------------------------------------------------- ##
    ## Sort Data based on time.id
    ## ---------------------------------------------------- ##
    oldTSCS <- cbind(time.id, subject.id, y, X, W)
    newTSCS <- oldTSCS[order(oldTSCS[,1]),]
    YT <- as.matrix(newTSCS[,3])
    XT <- as.matrix(newTSCS[,4:(4+K-1)])
    WT <- as.matrix(newTSCS[,(4+K):(4+K+Q-1)])

    nsubj  <- length(unique(subject.id))
    if (unique(subject.id)[1] != 1){
      stop("subject.id should start 1!")
    }

    ## subject.offset is the obs number from which a new subject unit starts
    subject.offset <- c(0, which(diff(sort(subject.id))==1)[-nsubj])
    ## col1: subj ID, col2: offset (C indexing), col3: #time periods in each subject
    nsubject.vec <- rep(NA, nsubj)
    for (i in 1:nsubj){
      nsubject.vec[i] <- sum(subject.id==unique(subject.id)[i])
    }
    subject.groupinfo <- cbind(unique(subject.id), subject.offset, nsubject.vec)

    ## time.groupinfo
    ## col1: time ID, col2: offset (C indexing), col3: # subjects in each time
    if(unique(time.id)[1] != 1){
      time.id <- time.id - unique(time.id)[1] + 1
      cat("time.id does not start from 1. So it is modified by subtracting the first unit of time.")
    }

    ntime <- max(nsubject.vec) ## maximum time length
    ntime.vec <- rep(NA, ntime)
    for (i in 1:ntime){
      ntime.vec[i] <- sum(time.id==unique(time.id)[i])
    }
    ## time.offset is the obs number from which a new time unit starts when we stack data by time.id
    time.offset <- c(0, which(diff(sort(time.id))==1)[-ntime])
    time.groupinfo <- cbind(unique(time.id), time.offset, ntime.vec)

    ## ---------------------------------------------------- ##
    ## Prior setting
    ## ---------------------------------------------------- ##
    ## mvn.prior <- MCMCpack:::form.mvn.prior(b0, B0, K)
    ## b0 <- mvn.prior[[1]]
    ## B0 <- mvn.prior[[2]]
    R0 <- as.matrix(R0)
    if(ncol(R0) != ncol(W)){
        stop("The dimension of R0 does not match the rank of W!")
    }
    ## B0inv   <-  solve(B0)
    R0inv   <-  solve(R0)

    ## prior inputs
    if (m > 0){
        P0 <- MCMCpack:::trans.mat.prior(m=m, n=ntime, a=a, b=b)
        ## initial values
        P  <- MCMCpack:::check.P(P.start, m, a=a, b=b)
    }else {
        P <- P0 <- matrix(1, 1, 1)
    }

    ## ---------------------------------------------------- ##
    ## Initialize.
    ## ---------------------------------------------------- ##
    ols    <-  lm(y ~ X-1)
    lambda <- rmvnorm(ns, rep(1, K));
    tau    <- rep(1, ns);
    alpha  <- rep(alpha[1], ns);

    if (!is.null(beta.start)) {
        beta <- matrix(beta.start, nrow = ns, ncol = length(beta.start))
    }else{
        cat("Initializing betas by SLOG\n")
        beta_slog <- SLOG(x = X, y = y, l = tau[1])
        beta <- matrix(beta_slog, nrow = ns, ncol = length(beta_slog))
    }

    if (is.na(sigma2.start[1])) {
        sig2  <- rep(summary(ols)$sigma, ns)^2
    } else{
        sig2 <- rep(sigma2start, ns)
    }

    ## ---------------------------------------------------- ##
    # initialize saving matrix
    ## ---------------------------------------------------- ##
    nstore <- mcmc/thin
    alphadraws  <- matrix(data=0, nstore, ns)
    taudraws    <- matrix(data=0, nstore, ns)
    betadraws   <- matrix(data=0, nstore, ns*K)
    lambdadraws <- matrix(data=0, nstore, ns*K)
    sigmadraws  <- matrix(data=0, nstore, ns)
    if(!fixed){
        Ddraws <- matrix(data=0, nstore, ns*Q*Q)
    }
    psdraws     <- matrix(data=0, ntime, ns)
    sdraws      <- matrix(data=0, nstore, ntime)
    Z.loglike.array <- matrix(data=0, nstore, NT)
    if (n.break > 0) Pmat <- matrix(NA, nstore, ns)

    known.alpha = FALSE
    XVX <- XVy <- ehat <- D <- Dinv <- br <- bi <- as.list(rep(NA, ns))
    for(j in 1:ns){
        XVX[[j]] <- matrix(0, K, K)
        XVy[[j]] <- matrix(0, K, 1)
        D[[j]] <- R0
        Dinv[[j]] <- R0inv
        br[[j]] <- matrix(NA, Q, N)
        bi[[j]] <- matrix(rnorm(Q), Q, 1)
     }

    ## make an array
    Yt_arr <- Xt_arr <- Wt_arr <- as.list(rep(NA, ntime))
    for (tt in 1:ntime){
        N.tt <- sum(time.id == tt)
        Yt_arr[[tt]] <- as.matrix(y[time.id == tt],  N.tt, 1)
        Xt_arr[[tt]] <- as.matrix(X[time.id == tt, ],  N.tt, K)
        Wt_arr[[tt]] <- as.matrix(W[time.id == tt, ],  N.tt, Q)
    }

    ## prepare initial stuff
    ## ps.store <- matrix(0, T, ns)
    start.time <- proc.time();
    ps.store   <- rep(0, ntime)
    totiter    <- mcmc+burn
    if (n.break > 0) {
        state <- sort(sample(1:ns, size=T, replace=TRUE, prob=(rep(1, ns))))
        ps <- matrix(1, T, ns)/ns
        ## cat("randomly chosen initial state = ", table(state), "\n")
    } else {
        state <- rep(1, T)
        ps <- matrix(1, T, 1)
    }

    if(verbose !=0){
        cat("----------------------------------------------------\n")
        cat("MCMC SparseChangeMixedPanel Sampler Starts! \n")
        cat("Initial state = ", table(state), "\n")
        ## cat("function called: ")
        ## print(call)
        ## cat("start.time: ", start.time, "\n")
        cat("----------------------------------------------------\n")
    }


    ## ---------------------------------------------------- ##
    ## MCMC sampler starts here!
    ## ---------------------------------------------------- ##
    random.perturb <- 0
    for(iter in 1:totiter){

        ## if( i%%verbose==0 ) cat("iteration ", i, "\n")
        if(iter == (burn+1) ) {
            ess.time = proc.time();
        }
        ## ---------------------------------------------------- ##
        ## Step 1: tau -- (no change)
        ## ---------------------------------------------------- ##
        ## for (j in 1:ns){
        ##     tau[j] = draw.tau(beta[j,], alpha[j], nu.shape, nu.rate)
        ##     ## cat("tau[",j, "] = ", tau[j], "\n");
        ## }

        tau <- draw_tau_cpp(beta, alpha, nu.shape, nu.rate, ns)

        ## ---------------------------------------------------- ##
        ## Step 2: bi (added)
        ## ---------------------------------------------------- ##
        SSE <- Nj <- rep(0, ns)
        for(j in 1:ns){
            XVX[[j]] <- matrix(0, K, K)
            XVy[[j]] <- matrix(0, K, 1)
        }

        if(fixed){
            for (j in 1:ns){
                ej <- state == j
                for (tt in which(ej)){
                    XVX[[j]] <- XVX[[j]] + crossprod(Xt_arr[[tt]])
                    XVy[[j]] <- XVy[[j]] + t(Xt_arr[[tt]])%*%Yt_arr[[tt]]
                }
            }
            for (i in 1:N){
                for (j in 1:ns){
                    ej <- as.numeric(is.element(time.id, which(state == j)) & subject.id==i)
                    nj <- sum(ej)
                    yj  <-  matrix(y[ej==1], nj, 1)
                    Xj  <-  matrix(X[ej==1,], nj, K)
                    ehatj <- yj - Xj%*%beta[j,]
                    ## Vj  <-  solve(sig2[j]*diag(nj))
                    ## XVX[[j]] <-  XVX[[j]] + crossprod(Xj)## t(Xj)%*%Xj
                    ## XVy[[j]] <-  XVy[[j]] + t(Xj)%*%yj
                    e = t(ehatj)%*%(ehatj);
                    SSE[j] = SSE[j] + e
                }
            }
        }else{
            for (i in 1:N){
                for (j in 1:ns){
                    ej <- as.numeric(is.element(time.id, which(state == j)) & subject.id==i)
                    nj <- sum(ej)
                    yj  <-  matrix(y[ej==1], nj, 1)
                    Xj  <-  matrix(X[ej==1,], nj, K)
                    Wj  <-  matrix(W[ej==1,], nj, Q)
                    ehatj <- yj - Xj%*%beta[j,]
                    Vj  <-   chol2inv(chol(sig2[j]*diag(nj) + Wj%*%D[[j]]%*%t(Wj)))
                    XVX[[j]] <-  XVX[[j]] + t(Xj)%*%Vj%*%Xj
                    XVy[[j]] <-  XVy[[j]] + t(Xj)%*%Vj%*%yj
                    V  <- chol2inv(chol(Dinv[[j]] + t(Wj)%*%Wj/sig2[j]))
                    U <- chol(V)
                    Mu <- V%*%(t(Wj)%*%ehatj)/sig2[j]
                    bi[[j]]      <- drop(Mu + t(U) %*% rnorm(Q)) ## as.vector(rmvnorm(1, post.bi.mean, post.bi.var))
                    br[[j]][,i]  <- bi[[j]]          # save all bi by subjectwise
                    ## wbr[j:(j+ni-1)]     <-  Wi%*%bi

                    ## FOR SIGMA
                    ## Nj[j] = YN[j] + yj.rows();
                    e = t(ehatj - Wj%*%bi[[j]])%*%(ehatj - Wj%*%bi[[j]]);
                    SSE[j] = SSE[j] + e
                    ## cat("sse: ", j, " = ", SSE[j], "\n")
                }
            }
        }
        ## cat("sse:  = ", SSE, "\n")

        ## bi_out <- draw_bi_cpp(y, X, W, D, Dinv, XVX, XVy, SSE, sig2, beta,
        ##                         state, time.id, subject.id, N, ns)
        ## XVy <- lapply(bi_out$XVy, function(x) matrix(x, K, 1))
        ## XVX <- lapply(bi_out$XVX, function(x) matrix(x, K, K, byrow = TRUE))
        ## SSE <- bi_out$SSE
        ## br  <- lapply(bi_out$br, function(x) matrix(x, Q, N, byrow = TRUE))
        ## cat("sse:  = ", SSE, "\n")

        ## ---------------------------------------------------- ##
        ## Step 3: sig2 (change!!)
        ## ---------------------------------------------------- ##
        for (j in 1:ns){
            Nj[[j]] <- sum(is.element(time.id, which(state == j)))
            shape = (c0 + Nj[[j]])*0.5;
            rate = (d0 + SSE[j])*0.5;
            sig2[j] = rinvgamma(1, shape, rate);
        }

        ## ---------------------------------------------------- ##
        ## Step 4: D (added)
        ## ---------------------------------------------------- ##
        if(!fixed){
            Nij <- Nj
            for (j in 1:ns){
                Nij[[j]] <- sum(is.element(1:max(time.id), which(state == j)))
                brbr    <- br[[j]]%*%t(br[[j]])
                Rinv    <- chol2inv(chol(R0inv + brbr))
                r       <- r0 + Nj[[j]] ## bi is repeated n.id times
                Dinv[[j]] <- rwish(r, Rinv)
                D[[j]]    <- chol2inv(chol(Dinv[[j]]))
            }
        }

        ## ---------------------------------------------------- ##
        ## Step 5: lambda (no change)
        ## ---------------------------------------------------- ##
        for (j in 1:ns){
            for(k in 1:K){
                lambda[j, k] =  2*retstable(0.5 * alpha[j], 1.0, beta[j, k]^2 / tau[j]^2, method="LD");
            }
        }

        ## ---------------------------------------------------- ##
        ## Step 6: beta (change)
        ## ---------------------------------------------------- ##
        beta <- draw_beta_svd_cpp(XVX, XVy, lambda, sig2, tau, ns, K)

        ## ---------------------------------------------------- ##
        ## Step 7: alpha (no change)
        ## ---------------------------------------------------- ##
        if (!known.alpha){
            for (j in 1:ns){
                if (!alpha.MH) {
                    alpha[j] <- draw.alpha.griddy(beta[j,], tau[j])
                } else {
                    alpha[j] <- draw.alpha(alpha[j], beta[j,], tau[j]);
                }
                ## cat("alpha[",j, "] = ", alpha[j], "\n");
            }
        }

        ## ---------------------------------------------------- ##
        ## Step 8: sampling S (change)
        ## ---------------------------------------------------- ##
        if (n.break > 0) {
            ## tau, alpha, lambda are all marginalized out in likelihood!
            ## state.out <- sparse.panel.state.sampler(m=m, T=T, N=N,  Yt_arr=Yt_arr, Xt_arr=Xt_arr,
            ##                                         Wt_arr=Wt_arr, beta=beta, bi=bi, sig2=sig2, D=D, P=P)
            ## JHP: we changed the mean of mvnormal distribution by adding br
            ## state.out <- sparse.panel.state.sampler2(m=m, T=T, N=N,  Yt_arr=Yt_arr, Xt_arr=Xt_arr,
            ##                                         beta=beta, br=br, sig2=sig2, D=D, P=P)
            if(fixed){
                state.out <- sparse_fixed_state_sampler_cpp(m = m, T = T, N = N, Yt_arr = Yt_arr,
                                                            Xt_arr = Xt_arr, beta = beta, sig2 = sig2, P = P)
            }else{
                state.out <- sparse_panel_state_sampler_cpp(m = m, T = T, N = N, Yt_arr = Yt_arr,
                                                            Xt_arr = Xt_arr, Wt_arr = Wt_arr, D = D,
                                                            beta = beta, sig2 = sig2, P = P)
            }
            state <- state.out$state
            ## report if random perturbation turns on only once.
            if(length(table(state)) < ns & random.perturb == 0){
                cat("\nThe number of sampled latent state is smaller than that of the designated state. This is due largely to model misfit. In order to do the model comparison, the latent state will be randomly perturbed by equal probabilities here. But we recommend users to rethink the model!\n")
                state <- sort(sample(1:ns, size=ntime, replace=TRUE, prob=apply(ps, 2, mean)))
                random.perturb <- random.perturb + 1
            }
            ps <- state.out$ps

            ## cat("table(state) = ", table(state),  "\n")
        }

        ## ---------------------------------------------------- ##
        ## Step 9: sampling P
        ## ---------------------------------------------------- ##
        if (n.break > 0) {
            switch_mat <- switchg(state, m = m)
            ## cat("switch = ", print(switch_mat) ,  "\n")
            for (j in 1:ns){
                switch1 <- P0[j,] + switch_mat[j,]
                pj      <- rdirichlet.cp(1, switch1)
                P[j,]   <- pj
            }
        }

        ## ---------------------------------------------------- ##
        ## report
        ## ---------------------------------------------------- ##
        if (verbose!= 0 &iter %% verbose == 0){
            cat("----------------------------------------------",'\n')
            cat("## iteration = ", iter, '\n')
            cat("----------------------------------------------",'\n')
            if (n.break > 0) {
                cat("sampled states: ", table(state), '\n')
                ## cat("Transition matrix: ", format(t(P), digits=2) , '\n')
                for(j in 1:ns){
                    cat("beta at state ", j, ": ", format(beta[j, ], digits=2), '\n')
                }
                if(!fixed){
                    for(j in 1:ns){
                        cat("D at state ", j, ": ", format(t(D[[j]]), digits=2) , '\n')
                    }
                }

            }else{
                cat("beta: ", format(t(beta), digits=2), '\n')
            }
            if(random.perturb>0){
                cat("random perturbation: ", random.perturb/iter, '\n')
            }
        }

        if (iter > burn && (iter %% thin == 0)) {
            alphadraws[(iter-burn)/thin,]  <- alpha
            betadraws[(iter-burn)/thin,]   <- t(beta)
            lambdadraws[(iter-burn)/thin,] <- t(lambda)
            sigmadraws[(iter-burn)/thin,]  <- sig2
            taudraws[(iter-burn)/thin,]    <- tau
            if(!fixed){
                Ddraws[(iter-burn)/thin,]      <- unlist(D)
            }
            if (n.break > 0) {
                Pmat[(iter-burn)/thin, ] <- diag(P)
                sdraws[(iter-burn)/thin,]      <- state
                ps.store <- ps.store + ps
            }

            if(Waic){
                marker <- 1
                for (i in 1:N){
                    for (j in 1:ns){
                        ej <- as.numeric(is.element(time.id, which(state == j)) & subject.id==i)
                        nj <- sum(ej)
                        if(nj > 0){
                            yj  <-  matrix(y[ej==1], nj, 1)
                            Xj  <-  matrix(X[ej==1,], nj, K)
                            
                            if(fixed){
                                mu.state  <-  Xj%*%beta[j,]
                            } else{
                                Wj  <-  matrix(W[ej==1,], nj, Q)
                                mu.state  <-  Xj%*%beta[j,] + Wj%*%bi[[j]]
                            }
                            ## cat("marker:(marker+nj-1) is ", marker:(marker+nj-1), "\n")
                            ## cat("dnorm is ", dnorm(yj, mean = mu.state, sd=sqrt(sig2[j]), log=TRUE), "\n")
                            Z.loglike.array[(iter-burn)/thin, marker:(marker+nj-1)] <-
                                dnorm(yj, mean = mu.state, sd=sqrt(sig2[j]), log=TRUE)
                            
                            ## log(dnorm(yi, Mu, sqrt(sig2)))
                        }
                        marker   <-  marker + nj
                        ## c(mu.state[t, state[t]])
                    }
                }
            }
        }
    }


    ## ---------------------------------------------------- ##
    ## Marginal Likelihood Estimation starts here!
    ## ---------------------------------------------------- ##
    ## ---------------------------------------------------- ##
    ## prepare
    ## ---------------------------------------------------- ##
    beta.st <- matrix(apply(betadraws, 2, mean), ns, K, byrow=TRUE)
    lambda.st <- matrix(apply(lambdadraws, 2, mean), ns, K, byrow=TRUE)
    sig2.st <- apply(sigmadraws, 2, mean)
    if(!fixed){
        D.st <- Dinv.st <- D
        D.st.raw <- apply(Ddraws, 2, mean)
        for(j in 1:ns){
            D.st[[j]] <- matrix(D.st.raw[(j-1)*Q*Q + 1:(Q*Q)], Q, Q)
            Dinv.st[[j]] <- chol2inv(chol(D.st[[j]]))
        }
    }
    
    alpha.st <- apply(alphadraws, 2, mean)
    tau.st <- apply(taudraws, 2, mean)
    if(n.break > 0){
        P.st <- apply(Pmat, 2, mean)
    }
    ## ---------------------------------------------------- ##
    ## Likelihood computation
    ## ---------------------------------------------------- ##
    marker <- 1
    resid <- loglike.t <- rep(NA, NT)
    bi.st <- bi
    for (i in 1:N){
        for (j in 1:ns){
            ej <- as.numeric(is.element(time.id, which(state == j)) & subject.id==i)
            nj <- sum(ej)
            if(nj > 0){
                yj  <-  matrix(y[ej==1], nj, 1)
                Xj  <-  matrix(X[ej==1,], nj, K)
                unscaled.yj  <-  matrix(unscaled.Y[ej==1], nj, 1)
                unscaled.Xj  <-  matrix(unscaled.X[ej==1,], nj, K)
                if(fixed){
                    mu.state.st  <-  Xj%*%beta.st[j,]
                    mu.state.unscaled  <-  unscaled.Xj%*%beta.st[j,]
                } else{
                    Wj  <-  matrix(W[ej==1,], nj, Q)
                    ehatj <- yj - Xj%*%beta.st[j,]
                    V  <- chol2inv(chol(Dinv[[j]] + t(Wj)%*%Wj/sig2[j]))
                    U <- chol(V)
                    Mu <- V%*%(t(Wj)%*%ehatj)/sig2[j]
                    bi.st[[j]]      <- drop(Mu + t(U) %*% rnorm(Q)) ## as.vector(rmvnorm(1, post.bi.mean, post.bi.var))
                    mu.state.st  <-  Xj%*%beta.st[j,] + Wj%*%bi.st[[j]]
                    ## unscaled.W must be considered here....later...JHP
                }
                ## Likelihood computation
                loglike.t[marker:(marker+nj-1)] <- dnorm(yj, mean = mu.state.st, sd=sqrt(sig2.st[j]), log=TRUE)
                resid[marker:(marker+nj-1)] <- yj - mu.state.st
            }
            marker   <-  marker + nj
        }
    }
    loglike <- sum(loglike.t)

    cat("\n---------------------------------------------- \n ")
    cat("Likelihood computation \n")
    cat("    loglike: ", as.numeric(loglike), "\n")
    cat("---------------------------------------------- \n ")
    
    if(marginal){
        ## holders
        density.sig2.holder <- density.nu.holder <- density.D.holder <- matrix(NA, mcmc, ns)

        ## ---------------------------------------------------- ##
        ## Marginal Step 1. density.nu
        ## ---------------------------------------------------- ##
        nu.st <- tau.st
        for(g in 1:mcmc){
            for(j in 1:ns){
                nu.st[j] <- tau.st[j]^(-alphadraws[g, j])
                shape.g <- nu.shape + K/alphadraws[g, j]
                rate.g <- nu.rate + sum(abs(betadraws[g,K*(j-1)+1:K])^alphadraws[g, j])
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
                Nj[[j]] <- sum(is.element(time.id, which(state == j)))
                shape = (c0 + Nj[[j]])*0.5;
                rate = (d0 + SSE[j])*0.5;
                sig2[j] = rinvgamma(1, shape, rate);
                density.sig2.holder[g, j] <- dinvgamma(sig2.st[j], shape, rate)
            }

            ## Fixed : tau
            ## tau <- draw_tau_cpp(beta, alpha, nu.shape, nu.rate, ns)

            ## Reduced Step 2: bi (added)
            SSE <- Nj <- rep(0, ns)
            for(j in 1:ns){
                XVX[[j]] <- matrix(0, K, K)
                XVy[[j]] <- matrix(0, K, 1)
            }

            if(fixed){
                for (j in 1:ns){
                    ej <- state == j
                    for (tt in which(ej)){
                        XVX[[j]] <- XVX[[j]] + crossprod(Xt_arr[[tt]])
                        XVy[[j]] <- XVy[[j]] + t(Xt_arr[[tt]])%*%Yt_arr[[tt]]
                    }
                }
                for (i in 1:N){
                    for (j in 1:ns){
                        ej <- as.numeric(is.element(time.id, which(state == j)) & subject.id==i)
                        nj <- sum(ej)
                        yj  <-  matrix(y[ej==1], nj, 1)
                        Xj  <-  matrix(X[ej==1,], nj, K)
                        ehatj <- yj - Xj%*%beta[j,]
                        e = t(ehatj)%*%(ehatj);
                        SSE[j] = SSE[j] + e
                    }
                }
            }else{
                for (i in 1:N){
                    for (j in 1:ns){
                        ej <- as.numeric(is.element(time.id, which(state == j)) & subject.id==i)
                        nj <- sum(ej)
                        yj  <-  matrix(y[ej==1], nj, 1)
                        Xj  <-  matrix(X[ej==1,], nj, K)
                        Wj  <-  matrix(W[ej==1,], nj, Q)
                        ehatj <- yj - Xj%*%beta[j,]
                        Vj  <-   chol2inv(chol(sig2[j]*diag(nj) + Wj%*%D[[j]]%*%t(Wj)))
                        XVX[[j]] <-  XVX[[j]] + t(Xj)%*%Vj%*%Xj
                        XVy[[j]] <-  XVy[[j]] + t(Xj)%*%Vj%*%yj
                        V  <- chol2inv(chol(Dinv[[j]] + t(Wj)%*%Wj/sig2[j]))
                        U <- chol(V)
                        Mu <- V%*%(t(Wj)%*%ehatj)/sig2[j]
                        bi[[j]]      <- drop(Mu + t(U) %*% rnorm(Q)) ## as.vector(rmvnorm(1, post.bi.mean, post.bi.var))
                        br[[j]][,i]  <- bi[[j]]          # save all bi by subjectwise
                        e = t(ehatj - Wj%*%bi[[j]])%*%(ehatj - Wj%*%bi[[j]]);
                        SSE[j] = SSE[j] + e
                    }
                }
            }
            ## Reduced Step 3: D (added)
            if(!fixed){
                Nij <- Nj
                for (j in 1:ns){
                    Nij[[j]] <- sum(is.element(1:max(time.id), which(state == j)))
                    brbr    <- br[[j]]%*%t(br[[j]])
                    Rinv    <- chol2inv(chol(R0inv + brbr))
                    r       <- r0 + Nij[[j]] ## bi is repeated n.id times
                    Dinv[[j]] <- rwish(r, Rinv)
                    D[[j]]    <- chol2inv(chol(Dinv[[j]]))
                }
            }


            ## Reduced Step 4: lambda (treating as latent)
            for (j in 1:ns){
                for(k in 1:K){
                    lambda[j, k] =  2*retstable(0.5 * alpha[j], 1.0, beta[j, k]^2 / tau.st[j]^2, method="LD");
                }
            }
            ## Reduced Step 5: beta
            beta <- draw_beta_svd_cpp(XVX, XVy, lambda, sig2, tau.st, ns, K)

            ## Reduced Step 6: alpha
            for (j in 1:ns){
                if (!alpha.MH) {
                    alpha[j] <- draw.alpha.griddy(beta[j,], tau.st[j])
                } else {
                    alpha[j] <- draw.alpha(alpha[j], beta[j,], tau.st[j]);
                }
            }
            ## Reduced Step 7: sampling S
            if (n.break > 0) {
                if(fixed){
                    state.out <- sparse_fixed_state_sampler_cpp(m = m, T = T, N = N, Yt_arr = Yt_arr,
                                                                Xt_arr = Xt_arr, beta = beta, sig2 = sig2, P = P)
                }else{
                    state.out <- sparse_panel_state_sampler_cpp(m = m, T = T, N = N, Yt_arr = Yt_arr,
                                                                Xt_arr = Xt_arr, Wt_arr = Wt_arr, D = D,
                                                                beta = beta, sig2 = sig2, P = P)
                }
                state <- state.out$state

                ## random perturbation
                if(length(table(state)) < ns & random.perturb == 0){
                    cat("\nThe number of sampled latent state is smaller than that of the designated state. This is due largely to model misfit. In order to do the model comparison, the latent state will be randomly perturbed by equal probabilities here. But we recommend users to rethink the model!\n")
                     state <- sort(sample(1:ns, size=ntime, replace=TRUE, prob=apply(ps, 2, mean)))
                }
            }

            ## Reduced Step 8: sampling P
            if(n.break > 0){
                switch  <-  switchg(state, m=m)
                for (j in 1:ns){
                    switch1 <-  P0[j,] + switch[j,]
                    pj      <-  rdirichlet.cp(1, switch1)
                    P[j,]  <-  pj
                }
            }
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
        if(FALSE) {
            for(g in 1:mcmc){
                ## Fixed sig2
                ## Fixed : tau

                ## Arrange panel data by a new state
                SSE <- Nj <- rep(0, ns)
                for(j in 1:ns){
                    Nj[[j]] <- sum(is.element(time.id, which(state == j)))
                    XVX[[j]] <- matrix(0, K, K)
                    XVy[[j]] <- matrix(0, K, 1)
                }

                ## Reduced Step 1: bi (added)
                if(fixed){
                    for (j in 1:ns){
                        ej <- state == j
                        for (tt in which(ej)){
                            XVX[[j]] <- XVX[[j]] + crossprod(Xt_arr[[tt]])
                            XVy[[j]] <- XVy[[j]] + t(Xt_arr[[tt]])%*%Yt_arr[[tt]]
                        }
                    }
                }else{
                    for (i in 1:N){
                        for (j in 1:ns){
                            ej <- as.numeric(is.element(time.id, which(state == j)) & subject.id==i)
                            nj <- sum(ej)
                            yj  <-  matrix(y[ej==1], nj, 1)
                            Xj  <-  matrix(X[ej==1,], nj, K)
                            Wj  <-  matrix(W[ej==1,], nj, Q)
                            ehatj <- yj - Xj%*%beta[j,]
                            Vj  <-   chol2inv(chol(sig2.st[j]*diag(nj) + Wj%*%D[[j]]%*%t(Wj)))
                            XVX[[j]] <-  XVX[[j]] + t(Xj)%*%Vj%*%Xj
                            XVy[[j]] <-  XVy[[j]] + t(Xj)%*%Vj%*%yj
                            V  <- chol2inv(chol(Dinv[[j]] + t(Wj)%*%Wj/sig2.st[j]))
                            U <- chol(V)
                            Mu <- V%*%(t(Wj)%*%ehatj)/sig2.st[j]
                            bi[[j]]      <- drop(Mu + t(U) %*% rnorm(Q)) ## as.vector(rmvnorm(1, post.bi.mean, post.bi.var))
                            br[[j]][,i]  <- bi[[j]]          # save all bi by subjectwise
                        }
                    }
                    ## Reduced Step 3: D (added)
                    Nij <- Nj
                    for (j in 1:ns){
                        Nij[[j]] <- sum(is.element(1:max(time.id), which(state == j)))
                        brbr    <- br[[j]]%*%t(br[[j]])
                        Rinv    <- chol2inv(chol(R0inv + brbr))
                        r       <- r0 + Nij[[j]] ## bi is repeated n.id times
                        Dinv[[j]] <- rwish(r, Rinv)
                        D[[j]]    <- chol2inv(chol(Dinv[[j]]))
                    }
                }

                ## Reduced Step 4: lambda (treating as latent)
                for (j in 1:ns){
                    for(k in 1:K){
                        lambda[j, k] =  2*retstable(0.5 * alpha[j], 1.0, beta[j, k]^2 / tau.st[j]^2, method="LD");
                    }
                }
                ## Reduced Step 5: beta
                for (j in 1:ns){
                    VInv = (XVX[[j]] + diag(lambda[j,] * sig2.st[j] / tau.st[j]^2, K));
                    V = chol2inv(chol(VInv));
                    U = chol(V) * sqrt(sig2.st[j]);
                    Mu = V %*% XVy[[j]];
                    beta[j,] = drop(Mu + t(U) %*% rnorm(K))
                    density.beta.holder[g, j]  <- log(dmvnorm(beta.st[j,], Mu, V))
                    ##     cat('beta [',j, ']', beta[j,], "\n")
                }
                ## beta <- draw_beta_cpp(XVX, XVy, lambda, sig2.st, tau.st, ns, K)
                ## beta <- draw_beta_cpp(XX, XY, lambda, sig2, tau.st, ns, K)

                ## Reduced Step 6: alpha
                for (j in 1:ns){
                    if (!alpha.MH) {
                        alpha[j] <- draw.alpha.griddy(beta[j,], tau.st[j])
                    } else {
                        alpha[j] <- draw.alpha(alpha[j], beta[j,], tau.st[j]);
                    }
                }
                ## Reduced Step 7: sampling S
                if (n.break > 0) {
                    if(fixed){
                        state.out <- sparse_fixed_state_sampler_cpp(m = m, T = T, N = N, Yt_arr = Yt_arr,
                                                                    Xt_arr = Xt_arr, beta = beta, sig2 = sig2.st, P = P)
                    }else{
                        state.out <- sparse_panel_state_sampler_cpp(m = m, T = T, N = N, Yt_arr = Yt_arr,
                                                                    Xt_arr = Xt_arr, Wt_arr = Wt_arr, D = D,
                                                                    beta = beta, sig2 = sig2.st, P = P)
                    }
                    state <- state.out$state
                    ## random sampling in case of missing states
                    if(length(table(state)) < ns & random.perturb == 0){
                        cat("\nThe number of sampled latent state is smaller than that of the designated state. This is due largely to model misfit. In order to do the model comparison, the latent state will be randomly perturbed by equal probabilities here. But we recommend users to rethink the model!\n")
                        state <- sort(sample(1:ns, size=ntime, replace=TRUE, prob=apply(ps, 2, mean)))
                    }

                }

                ## Reduced Step 8: sampling P
                if(n.break > 0){
                    switch  <-  switchg(state, m=m)
                    for (j in 1:ns){
                        switch1 <-  P0[j,] + switch[j,]
                        pj      <-  rdirichlet.cp(1, switch1)
                        P[j,]  <-  pj
                    }
                }
            }
            density.beta <- sum(log(apply(exp(density.beta.holder), 2, mean)))
            cat("\n---------------------------------------------- \n ")
            cat("Marignal Likelihood Computation Step 3 \n")
            cat("    density.beta: ", as.numeric(density.beta), "\n")
            cat("---------------------------------------------- \n ")
        }

        ## ---------------------------------------------------- ##
        ## Marginal Step 3: D.st|tau.st, sig2.st
        ## ---------------------------------------------------- ##
        if(n.break > 0 & !fixed){
            for(g in 1:mcmc){
                ## Fixed sig2
                ## Fixed tau
                ## Fixed beta

                ## Arrange panel data by a new state
                SSE <- Nj <- rep(0, ns)
                for(j in 1:ns){
                    Nj[[j]] <- sum(is.element(time.id, which(state == j)))
                    XVX[[j]] <- matrix(0, K, K)
                    XVy[[j]] <- matrix(0, K, 1)
                }

                ## Reduced Step 1: bi (added)
                for (i in 1:N){
                    for (j in 1:ns){
                        ej <- as.numeric(is.element(time.id, which(state == j)) & subject.id==i)
                        nj <- sum(ej)
                        yj  <-  matrix(y[ej==1], nj, 1)
                        Xj  <-  matrix(X[ej==1,], nj, K)
                        Wj  <-  matrix(W[ej==1,], nj, Q)
                        ehatj <- yj - Xj%*%beta.st[j,]
                        Vj  <-   chol2inv(chol(sig2.st[j]*diag(nj) + Wj%*%D[[j]]%*%t(Wj)))
                        XVX[[j]] <-  XVX[[j]] + t(Xj)%*%Vj%*%Xj
                        XVy[[j]] <-  XVy[[j]] + t(Xj)%*%Vj%*%yj
                        V  <- chol2inv(chol(Dinv[[j]] + t(Wj)%*%Wj/sig2.st[j]))
                        U <- chol(V)
                        Mu <- V%*%(t(Wj)%*%ehatj)/sig2.st[j]
                        bi[[j]]      <- drop(Mu + t(U) %*% rnorm(Q)) ## as.vector(rmvnorm(1, post.bi.mean, post.bi.var))
                        br[[j]][,i]  <- bi[[j]]          # save all bi by subjectwise
                    }
                }

                ## Reduced Step 2: D (added)
                Nij <- Nj
                for (j in 1:ns){
                    Nij[[j]] <- sum(is.element(1:max(time.id), which(state == j)))
                    brbr    <- br[[j]]%*%t(br[[j]])
                    Rinv    <- chol2inv(chol(R0inv + brbr))
                    r       <- r0 + Nij[[j]] ## bi is repeated n.id times
                    Dinv[[j]] <- rwish(r, Rinv)
                    D[[j]]    <- chol2inv(chol(Dinv[[j]]))
                    density.D.holder[g, j] <- log(dwish(Dinv[[j]], r, Rinv))
                }


                ## Reduced Step 3: lambda (treating as latent)
                for (j in 1:ns){
                    for(k in 1:K){
                        lambda[j, k] =  2*retstable(0.5 * alpha[j], 1.0, beta[j, k]^2 / tau.st[j]^2, method="LD");
                    }
                }
                ## Reduced Step 4: beta
                beta <- draw_beta_svd_cpp(XVX, XVy, lambda, sig2.st, tau.st, ns, K)
                ## beta <- draw_beta_cpp(XVX, XVy, lambda, sig2.st, tau.st, ns, K)
                ## beta <- draw_beta_cpp(XX, XY, lambda, sig2, tau.st, ns, K)

                ## Reduced Step 5: alpha
                for (j in 1:ns){
                    if (!alpha.MH) {
                        alpha[j] <- draw.alpha.griddy(beta[j,], tau.st[j])
                    } else {
                        alpha[j] <- draw.alpha(alpha[j], beta[j,], tau.st[j]);
                    }
                }
                ## Reduced Step 6: sampling S
                if(fixed){
                    state.out <- sparse_fixed_state_sampler_cpp(m = m, T = T, N = N, Yt_arr = Yt_arr,
                                                                Xt_arr = Xt_arr, beta = beta, sig2 = sig2.st, P = P)
                }else{
                    state.out <- sparse_panel_state_sampler_cpp(m = m, T = T, N = N, Yt_arr = Yt_arr,
                                                                Xt_arr = Xt_arr, Wt_arr = Wt_arr, D = D,
                                                                beta = beta, sig2 = sig2.st, P = P)
                }
                state <- state.out$state
                if(length(table(state)) < ns & random.perturb == 0){
                    cat("\nThe number of sampled latent state is smaller than that of the designated state. This is due largely to model misfit. In order to do the model comparison, the latent state will be randomly perturbed by equal probabilities here. But we recommend users to rethink the model!\n")
                    state <- sort(sample(1:ns, size=ntime, replace=TRUE, prob=apply(ps, 2, mean)))
                }


                ## Reduced Step 7: sampling P
                switch  <-  switchg(state, m=m)
                for (j in 1:ns){
                    switch1 <-  P0[j,] + switch[j,]
                    pj      <-  rdirichlet.cp(1, switch1)
                    P[j,]  <-  pj
                }

            }

            density.D <- sum(log(apply(exp(density.D.holder), 2, mean)))
            cat("\n---------------------------------------------- \n ")
            cat("Marignal Likelihood Computation Step 3 \n")
            cat("    density.D ", as.numeric(density.D), "\n")
            cat("---------------------------------------------- \n ")
        }

        ## ---------------------------------------------------- ##
        ## Marginal Step 4: P| tau.st, sig.st, D.st
        ## ---------------------------------------------------- ##
        if(n.break > 0){
            density.P.holder <- matrix(NA, mcmc, ns-1)
            for(g in 1:mcmc){
                ## Arrange panel data by a new state
                SSE <- Nj <- rep(0, ns)
                for(j in 1:ns){
                    Nj[[j]] <- sum(is.element(time.id, which(state == j)))
                    XVX[[j]] <- matrix(0, K, K)
                    XVy[[j]] <- matrix(0, K, 1)
                }
                if(fixed){
                    for (j in 1:ns){
                        ej <- state == j
                        for (tt in which(ej)){
                            XVX[[j]] <- XVX[[j]] + crossprod(Xt_arr[[tt]])
                            XVy[[j]] <- XVy[[j]] + t(Xt_arr[[tt]])%*%Yt_arr[[tt]]
                        }
                    }
                }else{
                    for (i in 1:N){
                        for (j in 1:ns){
                            ej <- as.numeric(is.element(time.id, which(state == j)) & subject.id==i)
                            nj <- sum(ej)
                            yj  <-  matrix(y[ej==1], nj, 1)
                            Xj  <-  matrix(X[ej==1,], nj, K)
                            Wj  <-  matrix(W[ej==1,], nj, Q)
                            ehatj <- yj - Xj%*%beta[j,]
                            Vj  <-   chol2inv(chol(sig2.st[j]*diag(nj) + Wj%*%D.st[[j]]%*%t(Wj)))
                            XVX[[j]] <-  XVX[[j]] + t(Xj)%*%Vj%*%Xj
                            XVy[[j]] <-  XVy[[j]] + t(Xj)%*%Vj%*%yj
                            ## V  <- chol2inv(chol(Dinv[[j]] + t(Wj)%*%Wj/sig2.st[j]))
                            ## U <- chol(V)
                            ## Mu <- V%*%(t(Wj)%*%ehatj)/sig2.st[j]
                            ## bi[[j]]      <- drop(Mu + t(U) %*% rnorm(Q)) ## as.vector(rmvnorm(1, post.bi.mean, post.bi.var))
                            ## br[[j]][,i]  <- bi[[j]]          # save all bi by subjectwise
                        }
                    }
                }

                ## Reduced Step 1: lambda
                for (j in 1:ns){
                    for(k in 1:K){
                        lambda[j, k] =  2*retstable(0.5 * alpha[j], 1.0, beta[j, k]^2 / tau.st[j]^2, method="LD");
                    }
                }
                ## Reduced Step 2: beta
                beta <- draw_beta_svd_cpp(XVX, XVy, lambda, sig2.st, tau.st, ns, K)

                ## Reduced Step 3: alpha
                for (j in 1:ns){
                    if (!alpha.MH) {
                        alpha[j] <- draw.alpha.griddy(beta[j,], tau.st[j])
                    } else {
                        alpha[j] <- draw.alpha(alpha[j], beta[j,], tau.st[j]);
                    }
                }
                ## Reduced Step 4: sampling S
                if(fixed){
                    state.out <- sparse_fixed_state_sampler_cpp(m = m, T = T, N = N, Yt_arr = Yt_arr,
                                                                Xt_arr = Xt_arr, beta = beta, sig2 = sig2.st, P = P)
                }else{
                    state.out <- sparse_panel_state_sampler_cpp(m = m, T = T, N = N, Yt_arr = Yt_arr,
                                                                Xt_arr = Xt_arr, Wt_arr = Wt_arr, D = D.st,
                                                                beta = beta, sig2 = sig2.st, P = P)
                }
                state <- state.out$state
                if(length(table(state)) < ns & random.perturb == 0){
                    cat("\nThe number of sampled latent state is smaller than that of the designated state. This is due largely to model misfit. In order to do the model comparison, the latent state will be randomly perturbed by equal probabilities here. But we recommend users to rethink the model!\n")
                    state <- sort(sample(1:ns, size=ntime, replace=TRUE, prob=apply(ps, 2, mean)))
                }


                ## Reduced Step  5: sampling P
                swit  <-  switchg(state, m=m)
                for (j in 1:ns){
                    swit1 <-  P0[j,] + swit[j,]
                    pj      <-  rdirichlet.cp(1, swit1)
                    P[j,]  <-  pj
                    if(j < ns){
                        shape1 <-  swit1[j]
                        shape2 <-  swit1[j + 1]
                        density.P.holder[g, j] <- dbeta(P.st[j], shape1, shape2)
                    }

                }
            }

            density.P <- log(prod(apply(density.P.holder, 2, mean)))
            cat("\n---------------------------------------------- \n ")
            cat("Marignal Likelihood Computation Step 5 \n")
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

        ## marginal prior density
        density.nu.prior <- density.sig2.prior <- density.D.prior <- rep(NA, ns)
        for (j in 1:ns){
            nu.st[[j]] <- tau.st[j]^(-alpha.st[j])
            density.nu.prior[j] <- log(dgamma(nu.st[[j]], nu.shape, rate = nu.rate))
            ## density.nu.prior[j] <- log(dgamma(nu.st[[j]], nu.shape, nu.rate))
            ## density.beta.prior[j] <- log(dmvnorm(beta.st[j, ], b0, B0))
            ## density.lambda.prior[j] <- log(dmvnorm(eV.st[[j]], eV0, iVV0))
            density.sig2.prior[j] <- log(dinvgamma(sig2.st[j], c0, d0))
            if(!fixed){
                density.D.prior[j] <- log(dwish(Dinv.st[[j]], r0, R0inv))
            }
        }

        ## ---------------------------------------------------- ##
        ## Log marginal likelihood addition
        ## ---------------------------------------------------- ##
         if(n.break > 0){
            if(fixed){
                logprior <- sum(density.nu.prior) + sum(density.sig2.prior) + sum(density.P.prior)
                logdenom <- (density.nu + density.sig2 + density.P)

            }else{
                logprior <- sum(density.nu.prior) + sum(density.sig2.prior) + sum(density.D.prior) + sum(density.P.prior)
                logdenom <- (density.nu + density.sig2 + density.P + density.D)
            }
           ## logmarglike <- (loglike + logprior) - logdenom;
        }else{
            if(fixed){
                logprior <- sum(density.nu.prior) + sum(density.sig2.prior)
                logdenom <- (density.nu + density.sig2)

            }else{
                logprior <- sum(density.nu.prior) + sum(density.sig2.prior) + sum(density.D.prior)
                logdenom <- (density.nu + density.sig2 + density.D)
            }
         }
        logmarglike <- (loglike + logprior) - logdenom

        cat("\n----------------------------------------------------\n")
        cat("    log marginal likelihood = (loglike + logprior) - (density.parameters) \n")
        cat("    log marginal likelihood : ", as.numeric(logmarglike), "\n")
        cat("    log likelihood : ", as.numeric(loglike), "\n")
        cat("    logprior : ", as.numeric(logprior), "\n")
        cat("    log posterior density : ", as.numeric(logdenom), "\n")
        cat("\n----------------------------------------------------\n")

    }## end of marginal

    end.time = proc.time();
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
    xnames <-  sapply(c(1:K), function(i){paste("beta", i, sep = "")})
    lnames <- sapply(c(1:K), function(i){paste("lambda", i, sep = "")})
    Dnames <-  sapply(c(1:(Q*Q)), function(i){paste("D", i, sep = "")})
    output <- NA
    if (m == 0){
        if (standardize) betadraws <- betadraws * ysd / Xsd

        if(fixed){
            output <- as.mcmc(cbind(betadraws, sigmadraws))
            colnames(output) <- c(xnames, "sigma2")
        }else{
            output <- as.mcmc(cbind(betadraws, sigmadraws, Ddraws))
            colnames(output) <- c(xnames, "sigma2", Dnames)
        }
    }else{
        if (standardize) {
            sidx <- rep(1:ns, each = ncol(X))
            xidx <- 1:ncol(betadraws)
            idx  <- split(xidx, sidx)
            C1   <- ysd / Xsd #sd(y) / apply(X, 2, sd)
            for (s in 1:ns) {
                betadraws[,idx[[s]]] <- t(apply(betadraws[,idx[[s]]], 1, function(x) x * C1))
            }
        }
        
        output1 <- coda::mcmc(data=betadraws, start=burn+1, end=burn + mcmc, thin=thin)
        output2 <- coda::mcmc(data=sigmadraws, start=burn+1, end=burn + mcmc, thin=thin)
        if(!fixed){
            output3 <- coda::mcmc(data=Ddraws, start=burn+1, end=burn + mcmc, thin=thin)
        }
        output4 <- coda::mcmc(data=lambdadraws, start=burn+1, end=burn + mcmc, thin=thin)
        colnames(output1)  <- sapply(c(1:ns),
                                     function(i){
                                         paste(xnames, "_regime", i, sep = "")
                                     })
        colnames(output2)  <- sapply(c(1:ns),
                                     function(i){
                                         paste("sigma2_regime", i, sep = "")
                                     })
        if(!fixed){
            colnames(output3)  <- sapply(c(1:ns),
                                         function(i){
                                             paste(Dnames, "_regime", i, sep = "")
                                         })
        }
        colnames(output4)  <- sapply(c(1:ns),
                                     function(i){
                                         paste(lnames, "_regime", i, sep = "")
                                     })

        if(!fixed){
            output    <- as.mcmc(cbind(output1, output2, output3))
        }else{
            output    <- as.mcmc(cbind(output1, output2))
        }
        ps.holder <- matrix(ps.store, ntime, ns)
        s.holder  <- matrix(sdraws, nstore, ntime)
    }
    ## attr(output, "X") <- X
    attr(output, "title") <- "SparseChangeMixedPanel Posterior Sample"
    ## attr(output, "call")   <- cl
    attr(output, "y")       <- y[1:ntime]
    attr(output, "X")       <- X[1:ntime, ]
    attr(output, "m")       <- m
    attr(output, "nsubj")   <- nsubj
    attr(output, "ntime")   <- ntime
    attr(output, "alpha")   <- coda::mcmc(data=alphadraws, start=burn+1, end=burn + mcmc, thin=thin)
    attr(output, "tau")     <- coda::mcmc(data=taudraws, start=burn+1, end=burn + mcmc, thin=thin)
    attr(output, "random.perturb")   <- random.perturb/totiter
    if (m > 0){
        attr(output, "s.store") <- s.holder
        prob.state <- cbind(sapply(1:ns, function(k){apply(s.holder == k, 2, mean)}))
        attr(output, "prob.state") <- prob.state
        attr(output, "lambda") <- output4
    }
    attr(output, "Waic.out") <- Waic.out
    attr(output, "loglike") <- loglike
    attr(output, "resid") <- resid
    if(marginal){
        attr(output, "logmarglike") <- logmarglike
    }
    class(output) <- c("mcmc", "BridgeChange")
    return(output)
}
