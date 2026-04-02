####################################
## no change
####################################
set.seed(1999)

## simulate data
K <- 100
n <- 100
X <- matrix(rnorm(n*K), n, K)
sig2 <- 4
beta.true <- matrix(rnorm(K), K, 1)*5
beta.true[sample(1:K, K/2, replace=FALSE)] <- 0
Y <- X%*%beta.true + rnorm(n, 0, sqrt(sig2))

#########################
## HMBB estimation
#########################
## Be aware that 100 runs are too short for analysis. 
mcmc = 100; burn = 100; verbose = 100; thin = 1;
formula <- Y ~ X
reg.cp0 <- BridgeChangeReg(formula=formula, n.break = 0, Waic = TRUE)
reg.cp1 <- BridgeChangeReg(formula=formula, n.break = 1, Waic = TRUE)
reg.cp2 <- BridgeChangeReg(formula=formula, n.break = 2, Waic = TRUE)

## model selection by WAIC
waic <- WaicCompare(list(reg.cp0, reg.cp1, reg.cp2), print = TRUE)
plotWaic(waic)

## obtain beta.hat
beta.est <- apply(reg.cp0[, grep("beta", colnames(reg.cp0))], 2, mean)

## plot
plot(beta.true, beta.est,  xlab="TRUE", ylab="EST", type = "n",
       xlim = range(beta.est), ylim = range(beta.true), asp = 1)
abline(a=0, b=1, col="red", lty = 3, lwd = 1.5)
points(beta.est, beta.true, col="darkblue")

####################################
## change-point case
####################################
set.seed(11199)
## simulate data
K <- 100
n <- 100
X <- matrix(rnorm(n*K), n, K)
sig2 <- 4
beta.true <- matrix(NA, 2, K)

beta.true[1,] <- matrix(rnorm(K, 1, 1), K, 1)*2
beta.true[2,] <- matrix(rnorm(K, -1, 1), K, 1)
# beta.true[1, sample(1:K, K/2, replace=FALSE)] <- 0
# beta.true[2, sample(1:K, K/2, replace=FALSE)] <- 0
mu1 <- X[1:(n/2), ]%*%beta.true[1,]
mu2 <- X[((n/2)+1):n, ]%*%beta.true[2,] 
Y  <- c(rnorm(n/2, mu1, sqrt(sig2)), rnorm(n/2, mu2, sqrt(sig2)))

## fit he model
formula <- Y ~ X
fit.cp0 <- BridgeChangeReg(formula=formula, mcmc=1000, burn=1000, n.break = 0, Waic = TRUE)
fit.cp1 <- BridgeChangeReg(formula=formula, mcmc=1000, burn=1000, n.break = 1, Waic = TRUE)
fit.cp2 <- BridgeChangeReg(formula=formula, mcmc=1000, burn=1000, n.break = 2, Waic = TRUE)

## model selection by WAIC
## WAIC prefers the no break model but the ground truth is the one break model
(waic <- WaicCompare(list(fit.cp0, fit.cp1, fit.cp2), print = TRUE))

## plot generated data
par(mfrow=c(1,3))
plot(1:length(beta.true[1,]), beta.true[1,],
    xlab="predictor", ylab="coefficients",
    ylim=range(beta.true), type='n', main="Parameter changes")
points(beta.true[1,], col="red", pch="1", cex=1)
points(beta.true[2,], col="blue", pch="2", cex=1)
plot(Y)
plotWaic(waic)

## true break at t = 50
par(mfrow=c(1, 2))
MCMCpack::plotState(fit.cp1, legend.control =c(60, 0.85), main="One break")
MCMCpack::plotState(fit.cp2, legend.control =c(60, 0.85), main="Two breaks")

## obtain beta.hat
beta.est <-  matrix(apply(fit.cp1[, grep("beta", colnames(fit.cp1))], 2, mean), 2, , byrow=TRUE)

## plot
par(mfrow=c(1, 2))
## regime 1
plot(beta.true[1,], beta.est[1,],  xlab="TRUE", ylab="EST", type = "n",main="Regime 1",
       xlim = range(beta.est), ylim = range(beta.true), asp = 1)
abline(a=0, b=1, col="red", lty = 3, lwd = 1.5)
points(beta.true[1,], beta.est[1,], col="darkblue")

## regime 2
plot(beta.true[2,], beta.est[2,],  xlab="TRUE", ylab="EST", type = "n", main="Regime 2",
       xlim = range(beta.est), ylim = range(beta.true), asp = 1)
abline(a=0, b=1, col="red", lty = 3, lwd = 1.5)
points(beta.true[2,], beta.est[2,], col="darkblue")

####################################
## dotplot over time
## time-varying movements of selected covariates
####################################
## all covariates
dotplotRegime(fit.cp1, hybrid=FALSE, location.bar=12, x.location="default",
              text.cex=0.8, main="Time-varying Movements of All Covariates")
