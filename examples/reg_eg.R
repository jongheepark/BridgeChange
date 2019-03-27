####################################
## no change
####################################

set.seed(1999)

## simulate data
K <- 100
n <- 100
X <- matrix(rnorm(n*K), n, K)
sig2 <- 0.2
beta.true <- matrix(rnorm(K), K, 1)*5
beta.true[sample(1:K, K/2, replace=FALSE)] <- 0
Y <- X%*%beta.true + rnorm(n, 0, sqrt(sig2))

## global parameter for estimation 
mcmc <- burn <- 100; thin <- 1; verbose <- 100;


## fit the model
out0 <- BridgeChangeReg(
    y = Y, X = X, n.break = 0,
    mcmc = mcmc, burn = burn,
    thin = thin, verbose = verbose
)

## obtain
beta.est <- coef(out0)

## plot
plot(beta.est, beta.true, ylab="TRUE", xlab="EST", type = "n",
       xlim = range(beta.est), ylim = range(beta.true), asp = 1)
abline(a=0, b=1, col="red", lty = 3, lwd = 1.5)
points(beta.est, beta.true, col="darkblue")

## summary of all results
## (show all coefficients and sigma estimates)
summary(out0)
plot(out0)


####################################
## change-point case
####################################

set.seed(1973);
## generate data
out <- BridgeChangeSim(
    ntime=100, predictor = 100,
    rho=0.0, constant.p =0,
    positive.jump=FALSE, varying.p = 0.5,
    break.point = 0.5, dgp.only=TRUE
)

## plot generated data
plot(1:length(out$true.beta[1,]), out$true.beta[1,],
    xlab="predictor", ylab="coefficients",
    ylim=range(out$true.beta), type='n')
points(out$true.beta[1,], col="red", pch="1", cex=1)
points(out$true.beta[2,], col="blue", pch="2", cex=1)

## fit he model
set.seed(11173);
out0 <- BridgeChangeReg(
    y = out$y.c, X = out$x.c, n.break = 0,
    scale.data=TRUE, intercept = TRUE,
    mcmc = mcmc, burn = burn, thin = thin, verbose = verbose,
    alpha.MH = TRUE, waic = TRUE
)

set.seed(11173);
out1 <- BridgeChangeReg(
    y = out$y.c, X = out$x.c, n.break = 1,
    scale.data = TRUE, intercept = TRUE,
    mcmc = mcmc, burn = burn, thin = thin, verbose = verbose,
    alpha.MH = TRUE, waic = TRUE
)


set.seed(11173);
out2 <- BridgeChangeReg(
    y=out$y.c, X=out$x.c, n.break = 2,
    scale.data = TRUE, intercept = TRUE,
    mcmc = mcmc, burn = burn, thin=thin, verbose=verbose,
    alpha.MH = TRUE, Waic = TRUE, marginal = FALSE
)

WaicCompare(list(out0, out1, out2))
# MarginalCompare(list(out0, out1))
