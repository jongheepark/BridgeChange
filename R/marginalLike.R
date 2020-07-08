<<<<<<< HEAD

marglike_BridgeChangeReg <- function() {


  ## ---------------------------------------------------- ##
  ## prepare
  ## ---------------------------------------------------- ##
  lambda.st <- matrix(apply(lambdadraws, 2, mean), ns, K, byrow=TRUE)
  sig2.st   <- apply(sigmadraws, 2, mean)
  alpha.st  <- apply(alphadraws, 2, mean)
  tau.st    <- apply(taudraws, 2, mean)
  state.st  <- round(apply(sdraws, 2, median))

  if (n.break > 0) {
    P.st <- apply(Pmat, 2, mean)
  }

  ## ---------------------------------------------------- ##
  ## Likelihood computation
  ## ---------------------------------------------------- ##
  loglike.t <- sapply(1:ntime, function(t){
    dnorm(ydm[t], mean = c(mu.st.state[t, state.st[t]]), sd = sqrt(sig2.st[state.st[t]]), log = TRUE)
  })

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
    

  return(list("loglike" = loglike, "logmarglike" = logmarglike))
}
=======
# 
# marglike_BridgeChangeReg <- function() {
# 
# 
#   ## ---------------------------------------------------- ##
#   ## prepare
#   ## ---------------------------------------------------- ##
#   lambda.st <- matrix(apply(lambdadraws, 2, mean), ns, K, byrow=TRUE)
#   sig2.st   <- apply(sigmadraws, 2, mean)
#   alpha.st  <- apply(alphadraws, 2, mean)
#   tau.st    <- apply(taudraws, 2, mean)
#   state.st  <- round(apply(sdraws, 2, median))
# 
#   if (n.break > 0) {
#     P.st <- apply(Pmat, 2, mean)
#   }
# 
#   ## ---------------------------------------------------- ##
#   ## Likelihood computation
#   ## ---------------------------------------------------- ##
#   loglike.t <- sapply(1:ntime, function(t){
#     dnorm(ydm[t], mean = c(mu.st.state[t, state.st[t]]), sd = sqrt(sig2.st[state.st[t]]), log = TRUE)
#   })
# 
#   loglike <- sum(loglike.t)
# 
#   cat("\n---------------------------------------------- \n ")
#   cat("Likelihood computation \n")
#   cat("    loglike: ", as.numeric(loglike), "\n")
#   cat("---------------------------------------------- \n ")
# 
#         ## ---------------------------------------------------- ##
#         ## holders
#         ## ---------------------------------------------------- ##
#         density.sig2.holder <- density.nu.holder <- density.beta.holder <- density.lambda.holder <- matrix(NA, mcmc, ns)
# 
# 
#         ## ---------------------------------------------------- ##
#         ## Marginal Step 1. density.nu
#         ## ---------------------------------------------------- ##
#         ## draw.tau <- function(beta, alpha, c, d)
#         ## {
#         ##   p = length(beta)
#         ##   nu = rgamma(1, c + p/alpha, rate=d + sum(abs(beta)^alpha))
#         ##   tau = nu^(-1/alpha)
#         ##   return(tau);
#         ## }
#         nu.st <- tau.st
#         for(g in 1:mcmc){
#             for(j in 1:ns){
#                 nu.st[j] <- tau.st[j]^(-alphadraws[g, j])
#                 shape.g  <- nu.shape + K/alphadraws[g, j]
#                 rate.g   <- nu.rate + sum(abs(betadraws[g,K*(j-1)+1:K])^alphadraws[g, j])
#                 density.nu.holder[g, j] <- dgamma(nu.st[j], shape.g, rate = rate.g)
#             }
#         }
#         density.nu <- log(prod(apply(density.nu.holder, 2, mean)))
#         cat("\n---------------------------------------------- \n ")
#         cat("Marignal Likelihood Computation Step 1 \n")
#         cat("    density.nu: ", as.numeric(density.nu), "\n")
#         cat("---------------------------------------------- \n ")
# 
#         ## ---------------------------------------------------- ##
#         ## Marginal Step 2: Sigma2|tau.st
#         ## ---------------------------------------------------- ##
#         for(g in 1:mcmc){
#             ## Reduced Step 1: sig2
#             for (j in 1:ns){
#                 ej  <-  as.numeric(state==j)
#                 nj  <-  sum(ej)
#                 yj  <-  matrix(ydm[ej==1], nj, 1)
#                 Xj  <-  matrix(X[ej==1,], nj, K)
#                 rss = sum( (as.matrix(yj)-Xj%*%as.matrix(beta[j,]))^2 )
#                 shape <- c0 + nj/2
#                 rate <- d0 + rss/2
#                 sig2[j] = rinvgamma(1, shape, scale=rate)
#                 density.sig2.holder[g, j] <- dinvgamma(sig2.st[j], shape, scale=rate)
#                 XX[[j]] <- t(Xj)%*%Xj
#                 XY[[j]] <- t(Xj)%*%yj
#             }
#             ## Fixed : tau
#             ## tau <- draw_tau_cpp(beta, alpha, nu.shape, nu.rate, ns)
# 
#             ## Reduced Step 2: lambda (treating as latent)
#             for (j in 1:ns){
#                 for(k in 1:K){
#                     lambda[j, k] =  2*retstable(0.5 * alpha[j], 1.0, beta[j, k]^2 / tau.st[j]^2, method="LD");
#                 }
#             }
#             ## Reduced Step 3: beta
#             beta <- draw_beta_svd_cpp(XX, XY, lambda, sig2, tau.st, ns, K)
#             ## beta <- draw_beta_cpp(XX, XY, lambda, sig2, tau.st, ns, K)
# 
#             ## Reduced Step 4: alpha
#             if (!known.alpha){
#                 if(alpha.limit){
#                     for (j in 1:ns){
#                         alpha[j] <- draw.alpha2(alpha[j], beta[j,], tau.st[j])
#                     }
#                 }else{
#                     ## if alpha is sampled from (0, 2]
#                     if (alpha.MH) {
#                         for (j in 1:ns){
#                             alpha[j] <- draw.alpha(alpha[j], beta[j,], tau.st[j]);
#                         }
#                         ## alpha[j] <- draw.alpha.randomwalk(alpha[j], beta[j,], tau[j], window=0.1)
#                     }else {
#                         for (j in 1:ns){
#                             alpha[j] <- draw.alpha.griddy(beta[j,], tau.st[j])
#                         }
#                     }
#                 }
#             }
# 
#             ## Reduced Step 5: sampling S
#             if(n.break > 0){
#                 state.out <- sparse_state_sampler_cpp(m, ydm, X, beta, sig2, P)
#                 state <- state.out$state
#                 ## if(length(table(state)) < ns){
#                 ##    state <- sort(sample(1:ns, size=ntime, replace=TRUE, prob=apply(ps, 2, mean)))
#                 ## }
# 
#             }
# 
#             ## Reduced Step 6: sampling P
#             if(n.break > 0){
#                 switch  <-  switchg(state, m=m)
#                 for (j in 1:ns){
#                     switch1 <-  A0[j,] + switch[j,]
#                     pj      <-  rdirichlet.cp(1, switch1)
#                     P[j,]  <-  pj
#                 }
#             }
# 
#             # ## demeaning y by regime
#             # if (n.break > 0 & demean == TRUE) {
#             #     ydm <- as.vector(as.vector(y) - tapply(y, state, mean)[state])
#             # }
# 
#         }
#         density.sig2 <- log(prod(apply(density.sig2.holder, 2, mean)))
#         cat("\n---------------------------------------------- \n ")
#         cat("Marignal Likelihood Computation Step 2 \n")
#         cat("    density.sig2: ", as.numeric(density.sig2), "\n")
#         cat("---------------------------------------------- \n ")
# 
#         ## ---------------------------------------------------- ##
#         ## Marginal Step 3: beta| tau.st, sig.st
#         ## Note that we do not evaluate the posterior ordinate for beta
#         ## as it is a latent variable with hyperparameter in Bayes Bridge model
#         ## ---------------------------------------------------- ##
#         ## for(g in 1:mcmc){
#         ##      if(n.break > 0){
#         ##         ## Fixed sig2
#         ##         for (j in 1:ns){
#         ##             ej  <-  as.numeric(state==j)
#         ##             nj  <-  sum(ej)
#         ##             yj  <-  matrix(ydm[ej==1], nj, 1)
#         ##             Xj  <-  matrix(X[ej==1,], nj, K)
#         ##             XX[[j]] <- t(Xj)%*%Xj
#         ##             XY[[j]] <- t(Xj)%*%yj
#         ##         }
#         ##     }
#         ## Fixed tau
#         ## tau <- draw_tau_cpp(beta, alpha, nu.shape, nu.rate, ns)
# 
#         ## Reduced Step 1: lambda
#         ##      for (j in 1:ns){
#         ##          for(k in 1:K){
#         ##              lambda[j, k] =  2*retstable(0.5 * alpha[j], 1.0, beta[j, k]^2 / tau.st[j]^2, method="LD");
#         ##          }
#         ##      }
#         ## Reduced Step 2: beta
#         ##      for (j in 1:ns){
#         ##          VInv = (XX[[j]] + diag(lambda[j,] * sig2.st[j] / tau.st[j]^2, K));
#         ##          V = chol2inv(chol(VInv));
#         ##          U = chol(V) * sqrt(sig2.st[j]);
#         ##          Mu = V %*% XY[[j]];
#         ##          beta[j,] = drop(Mu + t(U) %*% rnorm(K))
#         ##          density.beta.holder[g, j]  <- log(dmvnorm(beta.st[j,], Mu, V))
#         ##     cat('beta [',j, ']', beta[j,], "\n")
#         ##      }
#         ## Reduced Step 3: alpha
#         ##     for (j in 1:ns){
#         ##         if (!alpha.MH) {
#         ##              alpha[j] <- draw.alpha.griddy(beta[j,], tau.st[j])
#         ##         } else {
#         ##             alpha[j] <- draw.alpha(alpha[j], beta[j,], tau.st[j]);
#         ##         }
#         ##     }
#         ## Reduced Step 4: sampling S
#         ##     if(n.break > 0){
#         ##         state.out <- sparse_state_sampler_cpp(m, ydm, X, beta, sig2.st, P)
#         ##         state <- state.out$state
#         ##         if(length(table(state)) < ns){
#         ##             state <- sort(sample(1:ns, size=ntime, replace=TRUE, prob=apply(ps, 2, mean)))
#         ##         }
# 
#         ##         ## Reduced Step  5: sampling P
#         ##         switch  <-  switchg(state, m=m)
#         ##         for (j in 1:ns){
#         ##             switch1 <-  A0[j,] + switch[j,]
#         ##             pj      <-  rdirichlet.cp(1, switch1)
#         ##             P[j,]  <-  pj
#         ##         }
#         ##     }
#         ## }
#         ## density.beta <- sum(log(apply(exp(density.beta.holder), 2, mean)))
#         ## cat("\n---------------------------------------------- \n ")
#         ## cat("Marignal Likelihood Computation Step 3 \n")
#         ## cat("    density.beta: ", as.numeric(density.beta), "\n")
#         ## cat("---------------------------------------------- \n ")
# 
#         ## ---------------------------------------------------- ##
#         ## Marginal Step 3: P| tau.st, sig.st
#         ## ---------------------------------------------------- ##
#         if(n.break > 0){
#             density.P.holder <- matrix(NA, mcmc, ns-1)
#             for(g in 1:mcmc){
#                 if(n.break > 0){
#                     ## Fixed sig2
#                     for (j in 1:ns){
#                         ej  <-  as.numeric(state==j)
#                         nj  <-  sum(ej)
#                         yj  <-  matrix(ydm[ej==1], nj, 1)
#                         Xj  <-  matrix(X[ej==1,], nj, K)
#                         XX[[j]] <- t(Xj)%*%Xj
#                         XY[[j]] <- t(Xj)%*%yj
#                     }
#                 }
# 
#                 ## Reduced Step 1: lambda
#                 for (j in 1:ns){
#                     for(k in 1:K){
#                         lambda[j, k] =  2*retstable(0.5 * alpha[j], 1.0, beta[j, k]^2 / tau.st[j]^2, method="LD");
#                     }
#                 }
#                 ## Reduced Step 2: beta
#                 beta <- draw_beta_svd_cpp(XX, XY, lambda, sig2.st, tau.st, ns, K)
# 
#                 ## beta <- draw_beta_cpp(XX, XY, lambda, sig2.st, tau.st, ns, K)
# 
#                 ## Reduced Step 3: alpha
#                 if (!known.alpha){
#                     if(alpha.limit){
#                         for (j in 1:ns){
#                             alpha[j] <- draw.alpha2(alpha[j], beta[j,], tau.st[j])
#                         }
#                     }else{
#                         ## if alpha is sampled from (0, 2]
#                         if (alpha.MH) {
#                             for (j in 1:ns){
#                                 alpha[j] <- draw.alpha(alpha[j], beta[j,], tau.st[j]);
#                             }
#                             ## alpha[j] <- draw.alpha.randomwalk(alpha[j], beta[j,], tau[j], window=0.1)
#                         }else {
#                             for (j in 1:ns){
#                                 alpha[j] <- draw.alpha.griddy(beta[j,], tau.st[j])
#                             }
#                         }
#                     }
#                 }
#                 ## Reduced Step 4: sampling S
#                 state.out <- sparse_state_sampler_cpp(m, y, X, beta, sig2.st, P)
#                 state <- state.out$state
#                 if(length(table(state)) < ns){
#                     state <- sort(sample(1:ns, size=ntime, replace=TRUE, prob=apply(ps, 2, mean)))
#                 }
# 
# 
#                 ## Reduced Step  5: sampling P
#                 swit  <-  switchg(state, m=m)
#                 for (j in 1:ns){
#                     swit1 <-  A0[j,] + swit[j,]
#                     pj      <-  rdirichlet.cp(1, swit1)
#                     P[j,]  <-  pj
#                     if(j < ns){
#                         shape1 <-  swit1[j]
#                         shape2 <-  swit1[j + 1]
#                         density.P.holder[g, j] <- dbeta(P.st[j], shape1, shape2)
#                     }
#                 }
# 
#                 ## demeaning y by regime
#                 ## if (n.break > 0 & demean == TRUE) {
#                 ##     ydm <- as.vector(as.vector(y) - tapply(y, state, mean)[state])
#                 ## }
#             }
# 
#             density.P <- log(prod(apply(density.P.holder, 2, mean)))
#             cat("\n---------------------------------------------- \n ")
#             cat("Marignal Likelihood Computation Step 3 \n")
#             cat("    density.P: ", as.numeric(density.P), "\n")
#             cat("---------------------------------------------- \n ")
#         }
# 
#         ## ---------------------------------------------------- ##
#         ## prior density estimation
#         ## ---------------------------------------------------- ##
#         if(n.break > 0){
#             expected.duration <- round(ntime/(m + 1))
#             b <- 0.1
#             a <- b * expected.duration
#             density.P.prior <- rep(NA, ns-1)
#             for (j in 1:ns){
#                 if(j < ns){
#                     density.P.prior[j] <- log(dbeta(P.st[j], a, b)) ## p = 1
#                 }
#             }
#         }
# 
#         density.nu.prior <- density.sig2.prior <- rep(NA, ns)
#         for (j in 1:ns){
#             nu.st[[j]] <- tau.st[j]^(-alpha.st[j])
#             density.nu.prior[j] <- log(dgamma(nu.st[[j]], nu.shape, rate = nu.rate))
#             density.sig2.prior[j] <- log(dinvgamma(sig2.st[j], c0, d0))
#         }
# 
#         ## ---------------------------------------------------- ##
#         ## Log marginal likelihood addition
#         ## ---------------------------------------------------- ##
#         if(n.break > 0){
#             logprior <- sum(density.nu.prior) + sum(density.sig2.prior) + sum(density.P.prior)
#             logdenom <- (density.nu + density.sig2 + density.P)
#         }else{
#             logprior <- sum(density.nu.prior) + sum(density.sig2.prior)
#             logdenom <- (density.nu + density.sig2)
#         }
#         logmarglike <- (loglike + logprior) - logdenom
# 
# 
#         cat("\n----------------------------------------------------\n")
#         cat("    log marginal likelihood = (loglike + logprior) - (density.parameters) \n")
#         cat("    log marginal likelihood : ", as.numeric(logmarglike), "\n")
#         cat("    log likelihood : ", as.numeric(loglike), "\n")
#         cat("    logprior : ", as.numeric(logprior), "\n")
#         cat("    log.nu.prior : ", as.numeric(sum(density.nu.prior)), "\n")
#         cat("    log.sig2.prior : ", as.numeric(sum(density.sig2.prior)), "\n")
#         if(n.break > 0){cat("    log.P.prior : ", as.numeric(sum(density.P.prior)), "\n")}
#         cat("    log posterior density : ", as.numeric(logdenom), "\n")
#         cat("----------------------------------------------------\n")
#     }
# 
#   return(list("loglike" = loglike, "logmarglike" = logmarglike))
# }
>>>>>>> 9261927a984f545cb28f0202fdd59a6805067ac4
