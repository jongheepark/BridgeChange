#' Bridge Regression simulation
#'
#' Simulation code for univariate response change-point model with Bridge prior.
#'
#' @param N Number of cross-sectional units. If \code{N=1}, univariate time series data.
#' @param ntime Length of time series
#' @param predictor Number of predictor
#' @param rho correlation parameter 0 = no correlation.
#' @param time.series TRUE if dgp is generated from autocorrelated series. rho is used as an autocorrelation coefficient.
#' @param sign.change.tune tuning parameter for the size of parameter sign changes
#' @param sigma1 sigma 1
#' @param sigma2 sigma 2
#' @param train.ratio The proportion of training data. (0, 1).
#' @param fitted.mse If TRUE and n.break == 0, compute the MSE of fitted values against true responses. If FALSE, do the cross-validation test using training data.
#' @param corr.tune tuning parameter for sx
#' @param constant.p Proportion of constant parameters
#' @param varying.p Proportion of time-varying parameters
#' @param break.point Timing of a break between 0 and 1.
#' @param mcmc =100
#' @param burn =100
#' @param verbose =100 Verbose
#' @param thin Thinning
#' @param dgp.only If TRUE, only data are generated and estimation steps are skipped.
#' @return output
#'
#' @export
BridgeChangeSim <- function(ntime=500, predictor = 100, rho=0, time.series=FALSE, sign.change.tune = 2, sigma1=1, sigma2 = 2, train.ratio=0.5,
                            fitted.mse = TRUE, constant.p =0.1, varying.p = 0.2, break.point = 0.5, positive.jump=FALSE, n.break = 1, intercept=FALSE,
                            positive.jump.tune = 1, mcmc = 100, burn = 100, verbose = 100, thin = 1, N=1, known.alpha = FALSE,
                            dgp.only=FALSE){
    ## turn off warnings
    ## options(warn=-1)
    require(glmnet)
    m <- n.break; ns <- m + 1
    y <- rep(NA, ntime)
    cut <- ntime*break.point
    cons.predictor <- ceiling(constant.p * predictor)
    vary.predictor <- ceiling(varying.p * predictor)
    ## y <- apply(x[, 1:real_p], 1, sum) + rnorm(n)
    if(rho == 0){
        x <- matrix(rnorm(ntime*predictor), nrow=ntime, ncol=predictor)
    }
    if(rho != 0 & time.series){
        tmp.r <- matrix(rho, ntime, ntime)
        tmp.r <- tmp.r^abs(row(tmp.r)-col(tmp.r))
        x <- t(mvrnorm(predictor, rep(0, ntime), tmp.r))
    }
    if(rho != 0 & !time.series){
        require(mvtnorm)
        require(Matrix)

        ## no.corr.col <- sample(1:predictor, ceiling((1-rho)*predictor), prob=rep(1/predictor, predictor))
        R <- matrix(runif(predictor*predictor, max(0, rho-0.05), min(1, rho + 0.05)), ncol=predictor)
        ## R[, no.corr.col] <- 0
        R1 <- (R * lower.tri(R)) + t(R * lower.tri(R))
        ## diag(R1) <- 1
        diag(R1) <- 1 + abs(min(eigen(R1)$values))

        if(sum(eigen(R1)$values<0)>0){
            print(eigen(R1)$values)
            stop("Correlationa matrix is not positive semidefinite. ")

        }
        ## sigma <- riwish(P1+(1-rho)*P1, diag(P1)); image(sigma, main="Correlation Matrix")
        L = chol(R1)
        x <- t(t(L) %*% matrix(rnorm(ntime*predictor), nrow=predictor, ncol=ntime)) ## rmvnorm(ntime, rep(0, predictor), sigma=R1) ## matrix(rnorm(n*p), nrow=n, ncol=p)
        ## image(cor(x))
        ## x2 <- matrix(rnorm(T*P2), nrow=T, ncol=P2)
        ## x <- cbind(x1, x2)
    }
    true.beta <- matrix(NA, ns, predictor)
    if(n.break == 0){
        if(constant.p ==0){
            stop("constant.p must be larger than 0 when n.break = 0!")
        }else{
            permute <- sample(1:predictor, replace=FALSE)
            ## very small values for sparsity
            true.beta <- c(rnorm(cons.predictor, 0, 2), rep(0, predictor - cons.predictor))[permute]
            y <- x %*% true.beta + rnorm(ntime, 0, sigma1)
        }
    }else{
        if(cons.predictor > 1){
            if(positive.jump){ ## positive jump change
                y[1:cut] <- apply(x[1:cut, 1:cons.predictor], 1, sum) +
                    apply(sign.change.tune*x[1:cut, (cons.predictor+1):(cons.predictor+ vary.predictor)], 1, sum) + rnorm(cut, 0, sigma1)
                y[(cut+1):ntime] <- apply(x[(cut+1):ntime, 1:cons.predictor], 1, sum) +
                    apply(positive.jump.tune + sign.change.tune*x[(cut+1):ntime, (cons.predictor+1):(cons.predictor+ vary.predictor)], 1, sum) + rnorm(ntime - cut, 0, sigma2)
                true.beta[1, ] <- c(rep(1, cons.predictor), rep(sign.change.tune, vary.predictor), rep(0, predictor - vary.predictor - cons.predictor))
                true.beta[2, ] <- c(rep(1, cons.predictor), rep(positive.jump.tune + sign.change.tune, vary.predictor), rep(0, predictor - vary.predictor - cons.predictor))
            }else{ ## sign change only
                y[1:cut] <- apply(x[1:cut, 1:cons.predictor], 1, sum) +
                    apply(sign.change.tune*x[1:cut, (cons.predictor+1):(cons.predictor+ vary.predictor)], 1, sum) + rnorm(cut, 0, sigma1)
                y[(cut+1):ntime] <- apply(x[(cut+1):ntime, 1:cons.predictor], 1, sum) +
                    apply(-sign.change.tune*x[(cut+1):ntime, (cons.predictor+1):(cons.predictor+ vary.predictor)], 1, sum) + rnorm(ntime - cut, 0, sigma2)
                true.beta[1, ] <- c(rep(1, cons.predictor), rep(sign.change.tune, vary.predictor), rep(0, predictor - vary.predictor - cons.predictor))
                true.beta[2, ] <- c(rep(1, cons.predictor), rep(-sign.change.tune, vary.predictor), rep(0, predictor - vary.predictor - cons.predictor))
            }
        }else if(cons.predictor == 1){
            if(positive.jump){ ## positive jump change
                y[1:cut] <- sum(x[1:cut, cons.predictor]) +
                    apply(sign.change.tune*x[1:cut, (cons.predictor+1):(cons.predictor+ vary.predictor)], 1, sum) + rnorm(cut, 0, sigma1)
            y[(cut+1):ntime] <- sum(x[(cut+1):ntime, cons.predictor]) +
                apply(positive.jump.tune + sign.change.tune*x[(cut+1):ntime, (cons.predictor+1):(cons.predictor+ vary.predictor)], 1, sum) + rnorm(ntime - cut, 0, sigma2)
                true.beta[1, ] <- c(rep(1, cons.predictor), rep(sign.change.tune, vary.predictor), rep(0, predictor - vary.predictor - cons.predictor))
                true.beta[2, ] <- c(rep(1, cons.predictor), rep(positive.jump.tune + sign.change.tune, vary.predictor), rep(0, predictor - vary.predictor - cons.predictor))
            }else{ ## sign change only
                y[1:cut] <- sum(x[1:cut, cons.predictor]) +
                    apply(sign.change.tune*x[1:cut, (cons.predictor+1):(cons.predictor+ vary.predictor)], 1, sum) + rnorm(cut, 0, sigma1)
                y[(cut+1):ntime] <-  sum(x[(cut+1):ntime, cons.predictor])  +
                    apply(-sign.change.tune*x[(cut+1):ntime, (cons.predictor+1):(cons.predictor+ vary.predictor)], 1, sum) + rnorm(ntime - cut, 0, sigma2)
                true.beta[1, ] <- c(rep(1, cons.predictor), rep(sign.change.tune, vary.predictor), rep(0, predictor - vary.predictor - cons.predictor))
                true.beta[2, ] <- c(rep(1, cons.predictor), rep(-sign.change.tune, vary.predictor), rep(0, predictor - vary.predictor - cons.predictor))
            }
        }else{ ## cons.predictor == 0
            if(positive.jump){ ## positive jump in the slope ex) 2 -> 4
                y[1:cut] <- apply(sign.change.tune*x[1:cut, (cons.predictor+1):(cons.predictor+ vary.predictor)], 1, sum) + rnorm(cut, 0, sigma1)
                y[(cut+1):ntime] <- apply(positive.jump.tune + sign.change.tune*x[(cut+1):ntime, (cons.predictor+1):(cons.predictor+ vary.predictor)], 1, sum) + rnorm(ntime - cut, 0, sigma2)
                true.beta[1, ] <- c(rep(sign.change.tune, vary.predictor), rep(0, predictor - vary.predictor - cons.predictor))
                true.beta[2, ] <- c(rep(positive.jump.tune + sign.change.tune, vary.predictor), rep(0, predictor - vary.predictor - cons.predictor))
            }else{
                if( vary.predictor == 1){
                    y[1:cut] <- sum(sign.change.tune*x[1:cut, (cons.predictor+1):(cons.predictor+ vary.predictor)]) + rnorm(cut, 0, sigma1)
                    y[(cut+1):ntime] <- sum(-sign.change.tune*x[(cut+1):ntime, (cons.predictor+1):(cons.predictor+ vary.predictor)]) + rnorm(ntime - cut, 0, sigma2)
                    true.beta[1, ] <- c(rep(sign.change.tune, vary.predictor), rep(0, predictor - vary.predictor - cons.predictor))
                    true.beta[2, ] <- c(rep(-sign.change.tune, vary.predictor), rep(0, predictor - vary.predictor - cons.predictor))
                }else{
                    y[1:cut] <- apply(sign.change.tune*x[1:cut, (cons.predictor+1):(cons.predictor+ vary.predictor)], 1, sum) + rnorm(cut, 0, sigma1)
                    y[(cut+1):ntime] <- apply(-sign.change.tune*x[(cut+1):ntime, (cons.predictor+1):(cons.predictor+ vary.predictor)], 1, sum) + rnorm(ntime - cut, 0, sigma2)
                    true.beta[1, ] <- c(rep(sign.change.tune, vary.predictor), rep(0, predictor - vary.predictor - cons.predictor))
                    true.beta[2, ] <- c(rep(-sign.change.tune, vary.predictor), rep(0, predictor - vary.predictor - cons.predictor))
                }
            }
        }
    }
    ## shuffle x and y
    ## new.order <- sample(1:ntime, ntime)
    ## new.x <- x[new.order,]
    ## new.y <- y[new.order]

    ## Split data into train (.5) and test (.5) sets
    X.sd <- apply(x, 2, sd)
    y.sd <- sd(y)
    beta.true <-  true.beta*(X.sd/y.sd)

    train_rows <- sort(sample(1:ntime, train.ratio*ntime, prob=rep(1/ntime, ntime)))
    x.train <- x[train_rows, ]
    x.test <- x[-train_rows, ]

    y.train <- y[train_rows]
    y.test <- y[-train_rows]

    ## scale the data
    y.train.c <- scale(y.train)
    x.train.c <- apply(x.train, 2, scale)

    y.test.c <- scale(y.test)
    x.test.c <- apply(x.test, 2, scale)

    y.c <- scale(y)
    x.c <- apply(x, 2, scale)

    ## center the data
    ## y.train.c <- centerdata(matrix(y.train))
    ## x.train.c <- centerdata(x.train)

    ## y.test.c <- centerdata(matrix(y.test))
    ## x.test.c <- centerdata(x.test)

    ## y.c <- centerdata(matrix(y))
    ## x.c <- centerdata(x)

    ## plot
    # if(n.break > 0){
    #     mydata <- data.frame(
    #         x = c(apply(x.c[1:cut, ], 1, mean),apply(x.c[(cut+1):ntime, ], 1, mean)),
    #         y = y.c,
    #         Regime = as.factor(c(rep(1, length(1:cut)), rep(2, length((cut+1):ntime))))
    #     )
    #     print(ggplot(mydata, aes(x, y, fill = Regime)) +
    #           geom_smooth( aes(linetype = Regime, colour = Regime), method = "lm", ) +
    #           xlab("Mean of X") +
    #           ylab("Y") +
    #           geom_point( aes(shape = Regime, colour = Regime)))
    # }
    #
    # if(!dgp.only){
    #     if(n.break > 0){
    #         ## Fit the models
    #         ## Model 1 OLS
    #         fit.ols <- lm(y.train.c~x.train.c)
    #         yhat.ols <- predict(fit.ols, newx=x.test.c)
    #         mse.ols <- mean((y.test.c - yhat.ols)^2)
    #
    #         ## Model 2 Lasso
    #         fit.lasso <- cv.glmnet(x.train.c, y.train.c, type.measure="mse", alpha=1, standardize = TRUE, family="gaussian")
    #         yhat.lasso <- predict(fit.lasso, s=fit.lasso$lambda.1se, newx=x.test.c)
    #         mse.lasso <- mean((y.test.c - yhat.lasso)^2)
    #         ## plot(y.test.c, yhat.lasso); abline(a=0, b=1, col="brown")
    #
    #         ## Model 3 Ridge
    #         fit.ridge <- cv.glmnet(x.train.c, y.train.c, type.measure="mse", alpha=0, standardize = TRUE, family="gaussian")
    #         yhat.ridge <- predict(fit.ridge, s=fit.ridge$lambda.1se, newx=x.test.c)
    #         mse.ridge <- mean((y.test.c - yhat.ridge)^2)
    #         ## plot(y.test.c, yhat.ridge); abline(a=0, b=1, col="brown")
    #
    #         ## Model 4 Elastic
    #         fit.elastic <- cv.glmnet(x.train.c, y.train.c, type.measure="mse", alpha=0.5, standardize = TRUE, family="gaussian")
    #         yhat.elastic <- predict(fit.elastic, s=fit.elastic$lambda.1se, newx=x.test.c)
    #         mse.elastic <- mean((y.test.c - yhat.elastic)^2)
    #
    #         ## Model 5 Adaptive
    #         w3 <- 1/abs(matrix(coef(fit.ridge, s=fit.ridge$lambda.min)[, 1][2:(ncol(x)+1)] ))^1 ## Using gamma = 1
    #         w3[w3[,1] == Inf] <- 999999999 ## Replacing values estimated as Infinite for 999999999
    #         cv.lasso <- cv.glmnet(x=x.train.c, y=y.train.c, family='gaussian', alpha=1, penalty.factor=w3)
    #         yhat.adaptive <- predict(cv.lasso, s=cv.lasso$lambda.1se, newx=x.test.c)
    #         mse.adaptive <- mean((y.test.c - yhat.adaptive)^2)
    #
    #         ## Model 6 Fused lasso
    #         require(genlasso)
    #         fused.out = fusedlasso1d(y.train.c, X=x.train.c)
    #         beta.fused = coef(fused.out , lambda=1)$beta
    #         yfused <- predict(fused.out , Xnew=x.test.c, lambda=1)$fit
    #         mse.fused <- mean((y.test.c - yfused)^2)
    #
    #         ## Model 7 Lassoplus
    #         require(sparsereg)
    #         s1 <- sparsereg(y.train.c, X = x.train.c, EM=T)
    #         beta.lassoplus <- print(s1)[[1]][,1]
    #         yhat.lassoplus <- x.test.c%*%beta.lassoplus
    #         mse.lassoplus <- mean((y.test.c - yhat.lassoplus)^2)
    #
    #         ## model 8 SparseChange
    #         for (i in 0:n.break) {
    #             assign(paste("out", i, sep=""), SparseChangeReg(y=y.train.c, X=x.train.c, scale.data=TRUE, intercept = intercept,
    #                                                             mcmc=mcmc, beta.start = rnorm(predictor), demean=FALSE,
    #                                                             burn = burn, thin=thin, verbose=verbose, known.alpha = known.alpha,
    #                                                             n.break = i, Waic=TRUE, marginal=TRUE, alpha.MH = FALSE));
    #
    #             assign(paste("outMH", i, sep=""), SparseChangeReg(y=y.train.c, X=x.train.c, scale.data=TRUE, intercept=intercept,
    #                                                               mcmc=mcmc, beta.start = rnorm(predictor), demean=FALSE,
    #                                                               burn = burn, thin=thin, verbose=verbose, known.alpha = known.alpha,
    #                                                               n.break = i, Waic=TRUE, marginal=TRUE, alpha.MH = TRUE))
    #         }
    #         ## model.test <- WaicCompare(list(out0, out1, outMH0, outMH1))
    #         model.test <- -2*MarginalCompare(list(out0, out1, outMH0, outMH1))
    #         model.names <- c("nobreak_griddy", "onebreak_griddy", "nobreak_MH", "onebreak_MH")
    #         if(which.min(model.test) == 1 | which.min(model.test) == 3){
    #             beta.mat <- out0[, grep("beta", colnames(out0))]
    #             beta0 <- apply(beta.mat, 2, mean)
    #             mu0 <- x.test.c%*%beta0
    #             y.change <- c(mu0)
    #             mse.change0 <- mean((y.test.c - y.change)^2)
    #
    #             beta.mat <- outMH0[, grep("beta", colnames(out0))]
    #             beta0 <- apply(beta.mat, 2, mean)
    #             mu0 <- x.test.c%*%beta0
    #             y.change <- c(mu0)
    #             mse.changeMH0 <- mean((y.test.c - y.change)^2)
    #
    #             ## plot(attr(out0, "alpha"))
    #             mse.vec <- c(mse.ols, mse.lasso, mse.elastic, mse.ridge,
    #                          mse.adaptive, mse.fused, mse.lassoplus, mse.change0, mse.changeMH0)
    #             names(mse.vec) <- c("OLS", "Lasso", "Elastic", "Ridge",
    #                                 "Adaptive", "Fused", "LassoPlus", "SparseChangeGriddy", "SparseChangeMH")
    #
    #         }else{
    #             beta.mat <- out1[, grep("beta", colnames(out1))]
    #             beta1 <- apply(beta.mat[, grep("regime1", colnames(beta.mat))], 2, mean)
    #             beta2 <- apply(beta.mat[, grep("regime2", colnames(beta.mat))], 2, mean)
    #             state.vec <- apply(attr(out1, "s.store"), 2, median)
    #             mu1 <- x.test.c[state.vec==1, ]%*%beta1
    #             mu2 <- x.test.c[state.vec==2, ]%*%beta2
    #             y.change <- c(mu1, mu2)
    #             mse.change1 <- mean((y.test.c - y.change)^2)
    #             ## plot(attr(out1, "alpha"))
    #             beta.mat <- outMH1[, grep("beta", colnames(out1))]
    #             beta1 <- apply(beta.mat[, grep("regime1", colnames(beta.mat))], 2, mean)
    #             beta2 <- apply(beta.mat[, grep("regime2", colnames(beta.mat))], 2, mean)
    #             state.vec <- apply(attr(out1, "s.store"), 2, median)
    #             mu1 <- x.test.c[state.vec==1, ]%*%beta1
    #             mu2 <- x.test.c[state.vec==2, ]%*%beta2
    #             y.change <- c(mu1, mu2)
    #             mse.changeMH1 <- mean((y.test.c - y.change)^2)
    #             mse.vec <- c(mse.ols, mse.lasso, mse.elastic, mse.ridge,
    #                          mse.adaptive, mse.fused, mse.lassoplus, mse.change1, mse.changeMH1)
    #             names(mse.vec) <- c("OLS", "Lasso", "Elastic", "Ridge",
    #                                 "Adaptive", "Fused", "LassoPlus", "SparseChangeGriddy", "SparseChangeMH")
    #
    #         }
    #         cat("\n The chosen model is ", model.names[which.min(model.test)], "\n")
    #     } else{
    #         ## if n.break == 0 do the MSE test for all data
    #         ## Fit the models
    #         if(fitted.mse){
    #             ## Model 1 OLS
    #             ## fit.ols <- lm(y.c~x.c)
    #             ## yhat.ols <- predict(fit.ols, newx=x.c)
    #             ## mse.ols <- mean((y.c - yhat.ols)^2)
    #
    #             ## Model 2 Lasso
    #             fit.lasso <- cv.glmnet(x.c, y.c, type.measure="mse", alpha=1, standardize = TRUE, family="gaussian")
    #             yhat.lasso <- predict(fit.lasso, s=fit.lasso$lambda.1se, newx=x.c)
    #             mse.lasso <- mean((y.c - yhat.lasso)^2)
    #             ## plot(y.c, yhat.lasso); abline(a=0, b=1, col="brown")
    #
    #             ## Model 3 Ridge
    #             fit.ridge <- cv.glmnet(x.c, y.c, type.measure="mse", alpha=0, standardize = TRUE, family="gaussian")
    #             yhat.ridge <- predict(fit.ridge, s=fit.ridge$lambda.1se, newx=x.c)
    #             mse.ridge <- mean((y.c - yhat.ridge)^2)
    #             ## plot(y.c, yhat.ridge); abline(a=0, b=1, col="brown")
    #
    #             ## Model 4 Elastic
    #             fit.elastic <- cv.glmnet(x.c, y.c, type.measure="mse", alpha=0.5, standardize = TRUE, family="gaussian")
    #             yhat.elastic <- predict(fit.elastic, s=fit.elastic$lambda.1se, newx=x.c)
    #             mse.elastic <- mean((y.c - yhat.elastic)^2)
    #
    #             ## Model 5 Adaptive
    #             w3 <- 1/abs(matrix(coef(fit.ridge, s=fit.ridge$lambda.min)[, 1][2:(ncol(x)+1)] ))^1 ## Using gamma = 1
    #             w3[w3[,1] == Inf] <- 999999999 ## Replacing values estimated as Infinite for 999999999
    #             cv.lasso <- cv.glmnet(x=x.c, y=y.c, family='gaussian', alpha=1, penalty.factor=w3)
    #             yhat.adaptive <- predict(cv.lasso, s=cv.lasso$lambda.1se, newx=x.c)
    #             mse.adaptive <- mean((y.c - yhat.adaptive)^2)
    #
    #             ## Model 6 Fused lasso
    #             require(genlasso)
    #             fused.out = fusedlasso1d(y.c, X=x.c)
    #             beta.fused = coef(fused.out , lambda=1)$beta
    #             yfused <- predict(fused.out , Xnew=x.c, lambda=1)$fit
    #             mse.fused <- mean((y.c - yfused)^2)
    #
    #             ## Model 7 Lassoplus
    #             require(sparsereg)
    #             s1 <- sparsereg(y.c, X = x.c, EM=T)
    #             beta.lassoplus <- print(s1)[[1]][,1]
    #             yhat.lassoplus <- x.c%*%beta.lassoplus
    #             mse.lassoplus <- mean((y.c - yhat.lassoplus)^2)
    #
    #             model.names <- c("nobreak_griddy", "nobreak_MH")
    #             model <- c(2, 1)
    #             ## model 8 SparseChange
    #             out0 <- SparseChangeReg(y=y.c, X=x.c, scale.data=FALSE, intercept = intercept,
    #                                     mcmc=mcmc, beta.start = rnorm(predictor), demean=FALSE,
    #                                     burn = burn, thin=thin, verbose=verbose, known.alpha = known.alpha,
    #                                     n.break = 0, alpha.MH = FALSE);
    #
    #             outMH0 <- SparseChangeReg(y=y.c, X=x.c, scale.data=FALSE, intercept=intercept,
    #                                       mcmc=mcmc, beta.start = rnorm(predictor), demean=FALSE,
    #                                       burn = burn, thin=thin, verbose=verbose, known.alpha = known.alpha,
    #                                       n.break = 0, alpha.MH = TRUE)
    #
    #             beta.mat <- out0[, grep("beta", colnames(out0))]
    #             beta0 <- matrix(apply(beta.mat, 2, mean), predictor, 1)
    #             ## cat("dim x.c", "\n")
    #             mu0 <- x.c%*%beta0
    #             y.change0 <- c(mu0)
    #             mse.change0 <- mean((y.c - y.change0)^2)
    #
    #             beta.mat <- outMH0[, grep("beta", colnames(outMH0))]
    #             beta1 <- matrix(apply(beta.mat, 2, mean), predictor, 1)
    #             mu1 <- x.c%*%beta1
    #             y.change1 <- c(mu1)
    #             mse.changeMH0 <- mean((y.c - y.change1)^2)
    #         }else{
    #             ## Model 1 ols
    #             ## fit.ols <- lm(y.train.c~x.train.c)
    #             ## yhat.ols <- predict(fit.ols, newx=x.test.c)
    #             ## mse.ols <- mean((y.test.c - yhat.ols)^2)
    #
    #             ## Model 2 Lasso
    #             fit.lasso <- cv.glmnet(x.train.c, y.train.c, type.measure="mse", alpha=1, standardize = TRUE, family="gaussian")
    #             yhat.lasso <- predict(fit.lasso, s=fit.lasso$lambda.1se, newx=x.test.c)
    #             mse.lasso <- mean((y.test.c - yhat.lasso)^2)
    #             ## plot(y.test.c, yhat.lasso); abline(a=0, b=1, col="brown")
    #
    #             ## Model 3 Ridge
    #             fit.ridge <- cv.glmnet(x.train.c, y.train.c, type.measure="mse", alpha=0, standardize = TRUE, family="gaussian")
    #             yhat.ridge <- predict(fit.ridge, s=fit.ridge$lambda.1se, newx=x.test.c)
    #             mse.ridge <- mean((y.test.c - yhat.ridge)^2)
    #             ## plot(y.test.c, yhat.ridge); abline(a=0, b=1, col="brown")
    #
    #             ## Model 4 Elastic
    #             fit.elastic <- cv.glmnet(x.train.c, y.train.c, type.measure="mse", alpha=0.5, standardize = TRUE, family="gaussian")
    #             yhat.elastic <- predict(fit.elastic, s=fit.elastic$lambda.1se, newx=x.test.c)
    #             mse.elastic <- mean((y.test.c - yhat.elastic)^2)
    #
    #             ## Model 5 Adaptive
    #             w3 <- 1/abs(matrix(coef(fit.ridge, s=fit.ridge$lambda.min)[, 1][2:(ncol(x)+1)] ))^1 ## Using gamma = 1
    #             w3[w3[,1] == Inf] <- 999999999 ## Replacing values estimated as Infinite for 999999999
    #             cv.lasso <- cv.glmnet(x=x.train.c, y=y.train.c, family='gaussian', alpha=1, penalty.factor=w3)
    #             yhat.adaptive <- predict(cv.lasso, s=cv.lasso$lambda.1se, newx=x.test.c)
    #             mse.adaptive <- mean((y.test.c - yhat.adaptive)^2)
    #
    #             ## Model 6 Fused lasso
    #             require(genlasso)
    #             fused.out = fusedlasso1d(y.train.c, X=x.train.c)
    #             beta.fused = coef(fused.out , lambda=1)$beta
    #             yfused <- predict(fused.out , Xnew=x.test.c, lambda=1)$fit
    #             mse.fused <- mean((y.test.c - yfused)^2)
    #
    #             ## Model 7 Lassoplus
    #             require(sparsereg)
    #             s1 <- sparsereg(y.train.c, X = x.train.c, EM=T)
    #             beta.lassoplus <- print(s1)[[1]][,1]
    #             yhat.lassoplus <- x.test.c%*%beta.lassoplus
    #             mse.lassoplus <- mean((y.test.c - yhat.lassoplus)^2)
    #
    #             model.names <- c("nobreak_griddy", "nobreak_MH")
    #             model <- c(2, 1)
    #             ## model 8 SparseChange
    #             out0 <- BridgeChangeReg(y=y.train.c, X=x.train.c, scale.data=FALSE, intercept = intercept,
    #                                     mcmc=mcmc, beta.start = rnorm(predictor),
    #                                     burn = burn, thin=thin, verbose=verbose, known.alpha = known.alpha,
    #                                     n.break = 0, alpha.MH = FALSE);
    #
    #             outMH0 <- BridgeChangeReg(y=y.train.c, X=x.train.c, scale.data=FALSE, intercept=intercept,
    #                                       mcmc=mcmc, beta.start = rnorm(predictor),
    #                                       burn = burn, thin=thin, verbose=verbose, known.alpha = known.alpha,
    #                                       n.break = 0, alpha.MH = TRUE)
    #
    #             beta.mat <- out0[, grep("beta", colnames(out0))]
    #             beta0 <- matrix(apply(beta.mat, 2, mean), predictor, 1)
    #             ## cat("dim x.c", "\n")
    #             mu0 <- x.test.c%*%beta0
    #             y.change0 <- c(mu0)
    #             mse.change0 <- mean((y.test.c - y.change0)^2)
    #
    #             beta.mat <- outMH0[, grep("beta", colnames(outMH0))]
    #             beta1 <- matrix(apply(beta.mat, 2, mean), predictor, 1)
    #             mu1 <- x.test.c%*%beta1
    #             y.change1 <- c(mu1)
    #             mse.changeMH0 <- mean((y.test.c - y.change1)^2)
    #
    #         }
    #         ## plot(attr(out0, "alpha"))
    #         ## mse.vec <- c(mse.ols, mse.lasso, mse.elastic, mse.ridge,
    #         ##              mse.adaptive, mse.fused, mse.lassoplus, mse.change0, mse.changeMH0)
    #         ## names(mse.vec) <- c("OLS", "Lasso", "Elastic", "Ridge",
    #         ##                     "Adaptive", "Fused", "LassoPlus", "SparseChangeGriddy", "SparseChangeMH")
    #         mse.vec <- c(mse.lasso, mse.elastic, mse.ridge,
    #                      mse.adaptive, mse.fused, mse.lassoplus, mse.change0, mse.changeMH0)
    #         names(mse.vec) <- c("Lasso", "Elastic", "Ridge",
    #                             "Adaptive", "Fused", "LassoPlus", "BridgeGriddy", "BridgeMH")
    #
    #     }
    #     ## plot(mse.vec)
    #     library(ggplot2)
    #     theme_set(theme_bw())
    #
    #     ## Plot
    #     s <- data.frame(Model = names(mse.vec), MSE = mse.vec)
    #     sub.title <- paste0("T = ", ntime, " predictor = ", predictor, " rho = ", rho, " Varying predictor = ",
    #                         predictor*varying.p, " Constant predictor = ", constant.p*predictor )
    #     s$`Model name` <- rownames(s)
    #     s <- s[order(s$MSE, decreasing=TRUE), ]
    #     s$`Model name` <- factor(s$`Model name`, levels = s$`Model name`)
    #     gp <- ggplot(s, aes(x=`Model name`, y=MSE, label=MSE)) +
    #         geom_point(stat='identity', fill="black", size=6)  +
    #         geom_segment(aes(y = 0, x = `Model name`, yend = MSE, xend = `Model name`),
    #                          color = "black") +
    #         labs(title="MSE Comparison", subtitle=sub.title) +
    #         ylim(c(0, max(s$MSE))) + coord_flip()
    #     out <- list(mse=mse.vec, true.beta = true.beta,
    #                 y=y, x=x, y.c=y.c, x.c=x.c, y.sd = y.sd, X.sd = X.sd,
    #                 y.test.c=y.test.c, x.test.c=x.test.c,
    #                 y.train.c=y.train.c, x.train.c=x.train.c,
    #                 gp = gp, detected = model.names[which.min(model.test)])
    # }else{
        out <- list(y=y, x=x, y.c=y.c, x.c=x.c, y.sd = y.sd, X.sd = X.sd,
                    true.beta = beta.true,
                    y.test.c=y.test.c, x.test.c=x.test.c,
                    y.train.c=y.train.c, x.train.c=x.train.c)
    # }
    return(out)
    ## turn on warning
    ## options(warn=0)
}
