####################################
## no change
####################################
set.seed(1999)
mcmc = burn = 1000; thin=1; verbose=100;

## simulate data 
K <- 200
n <- 100
X <- matrix(rnorm(n*K), n, K)
sig2 <- .2
beta.hat <- matrix(rnorm(K), K, 1)*5
beta.hat[sample(1:K, K/2, replace=FALSE)] <- 0
Y <- X%*%beta.hat + rnorm(n, 0, sqrt(sig2))



out0 <- BridgeChangeReg(y = Y, X = X, 
                        scale.data = TRUE
                        mcmc = mcmc, burn = burn, thin = thin, verbose = verbose,
                        alpha.MH = TRUE, n.break = 0)
out11 <- BridgeChangeReg(y=Y, X=X, scale.data=FALSE,
                        mcmc=mcmc, intercept = TRUE,
                        burn = burn, thin=thin, verbose=verbose, known.alpha = 1,
                        alpha.MH=TRUE, n.break = 1)

X.sd <- apply(X, 2, sd)
y.sd <- sd(Y)
## normalized true.beta
beta.true <-  beta.hat*(X.sd/y.sd)

## estimate
beta.svd <- (apply(out0[, 1:200], 2, mean) + apply(X, 2, mean)) * apply(X, 2, sd)

## plot
plot(beta.svd, beta.true, ylab="TRUE", xlab="EST", type = "n",
       xlim = range(beta.svd), ylim = range(beta.svd))
points(beta.svd, beta.true, col="blue")
abline(a=0, b=1, col="red")

# ####################################
# ## one change with no correlation
# ####################################
# set.seed(1973);
# sim <- SparseChangeSim(ntime=100, predictor = 100, rho=0.0, constant.p = 0, sign.change.tune = 3, 
#                        positive.jump=FALSE, varying.p = 0.2, break.point = 0.5, dgp.only=TRUE)
# 
# plot(1:length(sim$true.beta[1,]), sim$true.beta[1,], xlab="predictor", ylab="coefficients", ylim=range(sim$true.beta), type='n')
# points(sim$true.beta[1,], col="red", pch="1", cex=1)
# points(sim$true.beta[2,], col="blue", pch="2", cex=1)
# 
# ## slope change
# set.seed(11173);
# out0 <- SparseChangeReg(y=sim$y.c, X=sim$x.c, scale.data=TRUE,
#                         mcmc=mcmc, demean=TRUE, intercept=FALSE,
#                         burn = burn, thin=thin, verbose=verbose,
#                         alpha.MH=TRUE, n.break = 0, Waic=TRUE, marginal = TRUE)
# 
# set.seed(11173);
# out1 <- SparseChangeReg(y=sim$y, X=sim$x, scale.data=TRUE,
#                         mcmc=100, demean=FALSE, intercept=FALSE,
#                         burn = 100, thin=1, verbose=verbose,
#                         alpha.MH=TRUE, n.break = 1, Waic=TRUE, marginal = TRUE)
# 
# WaicCompare(list(out0, out1))
# MarginalCompare(list(out0, out1))
# 
# SparseChange:::plotCoef(out1, scale.back = TRUE, true.beta = sim$true.beta)
# 
# 
# ## y=out$y.c; X=out$x.c; scale.data=TRUE;
# ## alpha.MH=TRUE; n.break = 0; Waic=TRUE; marginal = TRUE
# ## n.break = 1; scale.data = TRUE;
# ## mcmc = 100; burn = 100; verbose = 100; thin = 1;
# ## c0 = 0.1; d0 = 0.1; a = NULL; b = NULL; demean=FALSE; intercept=TRUE;
# ## beta.start = NULL; nu.shape=2.0; nu.rate=2.0; known.s = FALSE;
# ## alpha.start = 1; known.alpha = FALSE;
# 
# ## beta.st   <- matrix(apply(out1[, 1:400], 2, mean), ns, K, byrow=TRUE) ## matrix(apply(betadraws, 2, mean), ns, K, byrow=TRUE)
# ## lambda.st <- matrix(apply(attr(out1, "lambda"), 2, mean), ns, K, byrow=TRUE)
# ## sig2.st   <- apply(out1[, 401:402], 2, mean)
# ## alpha.st  <- apply(attr(out1, "alpha"), 2, mean)## apply(alphadraws, 2, mean)
# ## tau.st    <- apply(attr(out1, "tau"), 2, mean)## apply(taudraws, 2, mean)
# ## if(n.break > 0){
# ##     P.st <- apply(Pmat, 2, mean)
# ## }
# ## mu.st.state <- out$x.c %*% t(beta.st)
# 
# ####################################
# ## one change with small correlation
# ####################################
# set.seed(1973);
# mcmc = burn = 1000; thin=1; verbose=100;
# out <- SparseChangeSim(ntime=100, predictor = 120, rho=0.2, constant.p =0,
#                        positive.jump=FALSE, varying.p = 0.2, break.point = 0.5, dgp.only=TRUE)
# 
# plot(1:length(out$true.beta[1,]), out$true.beta[1,], xlab="predictor", ylab="coefficients", ylim=range(out$true.beta), type='n')
# points(out$true.beta[1,], col="red", pch="1", cex=1)
# points(out$true.beta[2,], col="blue", pch="2", cex=1)
# 
# set.seed(11173);
# out0 <- SparseChangeReg(y=out$y.c, X=out$x.c, scale.data=TRUE,
#                         mcmc=mcmc, demean=TRUE, intercept=FALSE,
#                         burn = burn, thin=thin, verbose=verbose,
#                         alpha.MH=TRUE, n.break = 0, Waic=TRUE, marginal = TRUE)
# 
# set.seed(11173);
# out1 <- SparseChangeReg(y=out$y.c, X=out$x.c, scale.data=TRUE,
#                         mcmc=mcmc, demean=TRUE, intercept=FALSE,
#                         burn = burn, thin=thin, verbose=verbose,  
#                         alpha.MH=TRUE, n.break = 1, Waic=TRUE, marginal = TRUE)
# 
# WaicCompare(list(out0, out1))
# MarginalCompare(list(out0, out1))
# 
# SparseChange:::plotCoef(out1, true.beta = out$true.beta)
# 
# ####################################
# ## one change with large correlation
# ####################################
# set.seed(1973);
# mcmc = burn = 1000; thin=1; verbose=100;
# out <- SparseChangeSim(ntime=100, predictor = 120, rho=0.5, constant.p =0,
#                        positive.jump=FALSE, varying.p = 0.2, break.point = 0.5, dgp.only=TRUE)
# 
# plot(1:length(out$true.beta[1,]), out$true.beta[1,], xlab="predictor", ylab="coefficients", ylim=range(out$true.beta), type='n')
# points(out$true.beta[1,], col="red", pch="1", cex=1)
# points(out$true.beta[2,], col="blue", pch="2", cex=1)
# 
# set.seed(11173);
# out0 <- SparseChangeReg(y=out$y.c, X=out$x.c, scale.data=TRUE,
#                         mcmc=mcmc, demean=TRUE, intercept=FALSE,
#                         burn = burn, thin=thin, verbose=verbose,
#                         alpha.MH=TRUE, n.break = 0, Waic=TRUE, marginal = TRUE)
# 
# set.seed(11173);
# out1 <- SparseChangeReg(y=out$y.c, X=out$x.c, scale.data=TRUE,
#                         mcmc=mcmc, demean=TRUE, intercept=FALSE,
#                         burn = burn, thin=thin, verbose=verbose,
#                         alpha.MH=TRUE, n.break = 1, Waic=TRUE, marginal = TRUE)
# 
# set.seed(11173);
# ## lasso-style constraint
# out11 <- SparseChangeReg(y=out$y.c, X=out$x.c, scale.data=TRUE,
#                         mcmc=mcmc, demean=TRUE, intercept=FALSE,
#                         burn = burn, thin=thin, verbose=verbose, known.alpha=1,
#                         alpha.MH=TRUE, n.break = 1, Waic=TRUE, marginal = TRUE)
# 
# WaicCompare(list(out0, out1, out11))
# MarginalCompare(list(out0, out1, out11))
# 
# SparseChange:::plotCoef(out11, true.beta = as.vector(out$true.beta))
# 
# ########################################################################
# ## one change with positive slope increases
# ########################################################################
# set.seed(1973);
# mcmc = burn = 1000; thin=1; verbose=100;
# out <- SparseChangeSim(ntime=100, predictor = 120, rho=0, constant.p =0,
#                        positive.jump=TRUE, varying.p = 0.2, break.point = 0.5, dgp.only=TRUE)
# 
# plot(1:length(out$true.beta[1,]), out$true.beta[1,], xlab="predictor", ylab="coefficients", ylim=range(out$true.beta), type='n')
# points(out$true.beta[1,], col="red", pch="1", cex=1)
# points(out$true.beta[2,], col="blue", pch="2", cex=1)
# 
# set.seed(11173);
# out0 <- SparseChangeReg(y=out$y.c, X=out$x.c, scale.data=TRUE,
#                         mcmc=mcmc, demean=TRUE, intercept=TRUE,
#                         burn = burn, thin=thin, verbose=verbose,
#                         alpha.MH=TRUE, n.break = 0, Waic=TRUE, marginal = TRUE)
# 
# set.seed(11173);
# out1 <- SparseChangeReg(y=out$y.c, X=out$x.c, scale.data=TRUE,
#                         mcmc=mcmc, demean=TRUE, intercept=TRUE,
#                         burn = burn, thin=thin, verbose=verbose,
#                         alpha.MH=TRUE, n.break = 1, Waic=TRUE, marginal = TRUE)
# 
# WaicCompare(list(out0, out1))
# MarginalCompare(list(out0, out1))
# 
# ## intercept added
# SparseChange:::plotCoef(out1, true.beta = c(0, 0, as.vector(out$true.beta)))
