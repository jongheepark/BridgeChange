######################################################################################################
## The idea is to develop a sparsity-induced prior-posterior model that fits data
## Bridge regression using mixture of normals representation.
## by JHP "Wed Oct 19 15:13:47 2016"
######################################################################################################

#' Hybrid Approach to Bridge Change Point Model with Fixed Effect
#'
#' @param fomula Inherited from \code{lm}. For example, \code{Y ~ X + Z}.
#' @param data Data.frame object.
#' @param index
#' String vector for unit and time index variables.
#' For example, \code{index = c("unit", "year")}.
#' @param model Model (\code{c("within","between", "pooling")}).
#' @param effect Effect (\code{c("individual", "time", "twoways")}).
#' @param standardize If TRUE, all covariates are standardized.
#' @param interaction If TRUE, all pairwise interactions are added.
#' @param n.break Number of breaks.
#' If \code{n.break = 0}, it simply runs fixed effect model with shrinkage prior on coefficients.
#' @param ols.weight If TRUE, OLS estimates are used for adpative lasso weight vector.
#' @param mcmc MCMC iteration.
#' @param burn Burn-in period.
#' @param verbose Verbose.
#' @param sparse.only If TRUE, skip MCMC and fit a sparse regressiong using adaptive lasso. This works only for n.break=0.
#' @param thin Thinning.
#' @param c0 Hyperparam
#' @param d0 = 0.1
#' @param nu.shape =2.0
#' @param nu.rate =2.0
#' @param alpha  = 1
#'
#'
#' @author Jong Hee Park, and Soichiro Yamauchi \email{syamauchi@princeton.edu}
#'
#' @export
#'
#'
## Adaptive for each regime
BridgeFixedPanelHybrid <- function(formula, data, index, model, effect,
                             standardize = TRUE, interaction = FALSE,
                             n.break = 1, ols.weight = FALSE, sparse.only = FALSE,
                             mcmc = 100, burn = 100, verbose = 100, thin = 1,
                             b0=0, B0=1, c0 = 0.1, d0 = 0.1, r0 =  1, R0 = 1,
                             nu.shape = 2.0, nu.rate = 2.0, alpha = 1,
                             Waic = FALSE, marginal = FALSE) {
    call <- match.call()
    a = NULL; b = NULL
    ## ---------------------------------------------------- ##
    ## use plm package here
    ## transform data int pdata.frame object
    ## ---------------------------------------------------- ##
    # if (standardize) {
    #     dat.sd <- apply(data[, !(colnames(data) %in% index)], 2, sd)
    #     data <- data.frame(cbind(scale(data[, !(colnames(data) %in% index)]), data[,index]))
    # }

    pdata    <- pdata.frame(data, index)
    pformula <- pFormula(formula)

    X <- plm:::model.matrix.pFormula(pformula, pdata, rhs = 1, model = model, effect = effect)
    y <- plm:::pmodel.response(pformula, pdata, model = model, effect = effect)

    plmX <- X
    plmy <- y

    m <- n.break
    W <- matrix(0, length(y), 1)

    if(interaction){
        x1.1 <- data.frame(X)
        var.names <- colnames(X)
        x1.2 <- matrix(t(apply(x1.1, 1, combn, 2, prod)), nrow = nrow(X))
        newX <- as.matrix(cbind(x1.1, x1.2))
        colnames(newX) <- c(var.names, combn(var.names, 2, paste, collapse="-"))
        X <- newX
    }
    unscaled.Y <- y
    unscaled.X <- X

    if (standardize) {
        ysd <- sd(y)
        Xsd <- apply(X, 2, sd)
        ## if Xsd is exactly 0, use 1.
        ## Xsd <- ifelse(Xsd == 0, 1, Xsd)
        dat.sd <- c(ysd, Xsd)
        X <- scale(X)
        y <- scale(as.vector(y))
    }


    subject.id <- as.numeric(as.factor(data[,index[1]]))
    time.id    <- as.numeric(as.factor(data[,index[2]]))

    ## ---------------------------------------------------- ##
    ## run change point model on demeaned data
    ## ---------------------------------------------------- ##
    if(sparse.only){
        output <- NA
    }else{
        output <- BridgeMixedPanel(subject.id = subject.id, time.id = time.id, y=as.vector(y), X=X, W=W,
                                   n.break = n.break, b0=b0, B0=B0, c0=c0, d0=d0, r0=r0, R0=R0,
                                   standardize = standardize, alpha.MH=TRUE,
                                   mcmc = mcmc, burn = burn, thin = thin, verbose=verbose,
                                   nu.shape = 2.0, nu.rate = 2.0, alpha = 1, Waic = Waic, marginal = marginal, fixed = TRUE,
                                   unscaled.Y = unscaled.Y, unscaled.X = unscaled.X)
    }
    ## data preparation for regime-wise regression
    ##  ninv.y <- 1/length(y)
    if(n.break > 0){
        state <- round(apply(attr(output, "s.store"), 2, mean))
        cat("estiamted states are ", table(state), "\n")
        beta.mean <- matrix(apply(output, 2, mean)[grepl("beta", colnames(output))], ncol(X), n.break+1) ## K by ns
        s2.mean <- matrix(apply(output, 2, mean)[grepl("sigma", colnames(output))], , n.break+1) ## 1 by ns
        yhat.state <- X%*%beta.mean
        ## yhat.state <- yhat.state - mean(yhat.state)
        state.indicator <- state[time.id]

        ## Following P. Richard HAHN and Carlos M. CARVALHO (Eq. 21)
        yhat <- sapply(1:length(y), function(tt){yhat.state[tt,state.indicator[tt]]})
        unique.time.index <- sort(unique(attr(plmX,"index")[,2]))
        n.state <- length(unique(state))
        raw.y.list <- y.list <- x.list <- as.list(rep(NA, n.state))
        ## raw.y.list.0 <- y.list.0 <- x.list.0 <- as.list(rep(NA, n.state))
        for(i in 1:n.state){
            ## x.list[[i]] <- X[ is.element(attr(plmX,"index")[,2], unique.time.index[state==i]), ]
            ## augmenting data frame with cluster-mean centered variables
            dat.x <- data.frame(X[ is.element(attr(plmX,"index")[,2], unique.time.index[state==i]), ])
            dat.y <- data.frame(yhat[ is.element(attr(plmX,"index")[,2], unique.time.index[state==i])])
            raw.y <- data.frame(y[ is.element(attr(plmX,"index")[,2], unique.time.index[state==i])])
            dat.id <- subject.id[is.element(attr(plmX,"index")[,2], unique.time.index[state==i])]

            ## group centering
            x.list[[i]] <- as.matrix(group.center(dat.x, dat.id))
            y.list[[i]] <- as.matrix(group.center(dat.y, dat.id))
            raw.y.list[[i]] <- as.matrix(group.center(raw.y, dat.id))
            ## x.list.0[[i]] <- as.matrix(dat.x)
            ## y.list.0[[i]] <- as.matrix(dat.y)
            ## raw.y.list.0[[i]] <- as.matrix(raw.y)
        }
        ## Variable selection using adaptive lasso
        if(ols.weight){
            ## Using OLS beta as weight vector
            hybrid.dss <- sapply(1:n.state, function(i){adaptive.lasso.olsweight(y.list[[i]], x.list[[i]])})
            rownames(hybrid.dss) <- colnames(X)
            colnames(hybrid.dss) <- paste0("Regime", 1:n.state)

            hybrid.cp <- sapply(1:n.state, function(i){adaptive.lasso.olsweight(raw.y.list[[i]], x.list[[i]])})
            rownames(hybrid.cp) <- colnames(X)
            colnames(hybrid.cp) <- paste0("Regime", 1:n.state)
        }else{
            hybrid.dss <- sapply(1:n.state, function(i){adaptive.lasso(y.list[[i]], x.list[[i]], beta.hat = beta.mean[,i])})
            rownames(hybrid.dss) <- colnames(X)
            colnames(hybrid.dss) <- paste0("Regime", 1:n.state)

            hybrid.cp <- sapply(1:n.state, function(i){adaptive.lasso(raw.y.list[[i]], x.list[[i]], beta.hat = beta.mean[,i])})
            rownames(hybrid.cp) <- colnames(X)
            colnames(hybrid.cp) <- paste0("Regime", 1:n.state)
        }
        ## hybrid.dss.0 <- sapply(1:n.state, function(i){adaptive.lasso(y.list.0[[i]], x.list.0[[i]])})
        ## rownames(hybrid.dss.0) <- colnames(X)
        ## colnames(hybrid.dss.0) <- paste0("Regime", 1:n.state)

        ## hybrid.cp.0 <- sapply(1:n.state, function(i){adaptive.lasso(raw.y.list.0[[i]], x.list.0[[i]])})
        ## rownames(hybrid.cp.0) <- colnames(X)
        ## colnames(hybrid.cp.0) <- paste0("Regime", 1:n.state)

        R2.DSS <- R2.CP <- R2.DSS.jhp <- R2.CP.jhp <- rep(NA, n.state)
        for(i in 1:n.state){
            R2.DSS[i] <- r.square.sparse(y = raw.y.list[[i]], yhat = y.list[[i]], x = x.list[[i]],
                                      beta.sparse = hybrid.dss[,i], s2.mean[i])
            R2.CP[i] <- r.square.sparse(y = raw.y.list[[i]], yhat = y.list[[i]], x = x.list[[i]],
                                      beta.sparse = hybrid.cp[,i], s2.mean[i])
            R2.DSS.jhp[i] <- r.square.sparse.jhp(y = raw.y.list[[i]], yhat = y.list[[i]], x = x.list[[i]],
                                      beta.sparse = hybrid.dss[,i])
            R2.CP.jhp[i] <- r.square.sparse.jhp(y = raw.y.list[[i]], yhat = y.list[[i]], x = x.list[[i]],
                                      beta.sparse = hybrid.cp[,i])
            ## cat("R-squared (DSS) : ", R2.DSS[i], "and R-squared (CP) ", R2.CP[i], "\n")
            cat("Pseudo R-squared (DSS) : ", R2.DSS.jhp[i], "and Pseudo R-squared (CP) ", R2.CP.jhp[i], "\n")
        }
        ## summary of pseudo R-squared
        yhat.dss <- X%*%hybrid.dss
        yhat.cp <- X%*%hybrid.cp
        yhat.sparse.dss <- sapply(1:length(y), function(tt){yhat.dss[tt, state.indicator[tt]]})
        yhat.sparse.cp <- sapply(1:length(y), function(tt){yhat.cp[tt, state.indicator[tt]]})
        R2.DSS.jhp.total <- r.square.sparse.total(y = y, yhat = yhat, yhat.sparse=yhat.sparse.dss)
        R2.CP.jhp.total <- r.square.sparse.total(y = y, yhat = yhat, yhat.sparse=yhat.sparse.cp)

        cat("Total :---------------------------------------------------------------------- \n Pseudo R-squared (DSS) : ", R2.DSS.jhp.total, "and Pseudo R-squared (CP) ", R2.CP.jhp.total, "\n")


    }else{
        if(sparse.only){
            dat.id <- subject.id
            ## group centering
            dat.x <- as.matrix(group.center(X, dat.id))
            raw.y <- as.matrix(group.center(matrix(y), dat.id))
            hybrid.cp <- matrix(adaptive.lasso.olsweight(raw.y, dat.x), ncol(X), 1)
            rownames(hybrid.cp) <- colnames(X)
            hybrid.dss <- NA
        }else{
            ## Following P. Richard HAHN and Carlos M. CARVALHO (Eq. 21)
            beta.mean <- matrix(apply(output, 2, mean)[grepl("beta", colnames(output))], ncol(X),1)
            s2.mean <- matrix(apply(output, 2, mean)[grepl("sigma", colnames(output))], , n.break+1) ## 1 by ns
            yhat <- X%*%beta.mean
            ## yhat <- yhat - mean(yhat)
            dat.id <- subject.id
            ## cat("y and yhat correlation is ", cor(yhat, y), "\n")

            ## group centering
            dat.x <- as.matrix(group.center(X, dat.id))
            dat.y <- as.matrix(group.center(matrix(yhat), dat.id))

            ## Variable selection using adaptive lasso
            if(ols.weight){
                hybrid.dss <- matrix(adaptive.lasso.olsweight(dat.y, dat.x), ncol(X), 1)
                rownames(hybrid.dss) <- colnames(X)

                raw.y <- as.matrix(group.center(matrix(y), dat.id))
                hybrid.cp <- matrix(adaptive.lasso.olsweight(raw.y, dat.x), ncol(X), 1)
                rownames(hybrid.cp) <- colnames(X)
            }else{
                hybrid.dss <- matrix(adaptive.lasso(dat.y, dat.x, beta.mean), ncol(X), 1)
                rownames(hybrid.dss) <- colnames(X)

                raw.y <- as.matrix(group.center(matrix(y), dat.id))
                hybrid.cp <- matrix(adaptive.lasso(raw.y, dat.x, beta.mean), ncol(X), 1)
                rownames(hybrid.cp) <- colnames(X)
            }
            i = 1
            R2.DSS <- r.square.sparse(y = raw.y, yhat = yhat, x = X,
                                      beta.sparse = hybrid.dss, s2.mean)
            R2.CP <- r.square.sparse(y = raw.y, yhat = yhat, x = X,
                                     beta.sparse = hybrid.cp, s2.mean)
            R2.DSS.jhp <- r.square.sparse.jhp(y = raw.y, yhat = yhat, x = X,
                                              beta.sparse = hybrid.dss)
            R2.CP.jhp <- r.square.sparse.jhp(y = raw.y, yhat = yhat, x = X,
                                             beta.sparse = hybrid.cp)
            cat("R-squared (DSS) : ", R2.DSS, "and R-squared (CP) ", R2.CP, "\n")

        }
    }
    attr(output, "hybrid") <- hybrid.dss
    attr(output, "hybrid.raw") <- hybrid.cp
    if(sparse.only){
    }else{
        attr(output, "R2") <- list("R2.DSS" = R2.DSS, "R2.CP" = R2.CP)
        attr(output, "R2.jhp") <- list("R2.DSS" = R2.DSS.jhp, "R2.CP" = R2.CP.jhp)
        if(n.break>0){
            attr(output, "R2.total") <- list("R2.DSS" = R2.DSS.jhp.total, "R2.CP" = R2.CP.jhp.total)
        }else{
            attr(output, "R2.total") <- list("R2.DSS" = R2.DSS.jhp, "R2.CP" = R2.CP.jhp)
        }
    }
    attr(output, "title")  <- "SparseChangeFixedPanel Posterior Sample"
    attr(output, "m")      <- n.break
    if(standardize) attr(output, "dat.sd") <- dat.sd
    return(output)
}

######################################################################################################
## The idea is to develop a sparsity-induced prior-posterior model that fits data
## Bridge regression using mixture of normals representation.
## by JHP "Wed Oct 19 15:13:47 2016"
######################################################################################################

#' Hybrid Approach to Bridge Change Point Model with Fixed Effect
#'
#' @param fomula Inherited from \code{lm}. For example, \code{Y ~ X + Z}.
#' @param data Data.frame object.
#' @param index
#' String vector for unit and time index variables.
#' For example, \code{index = c("unit", "year")}.
#' @param model Model (\code{c("within","between", "pooling")}).
#' @param effect Effect (\code{c("individual", "time", "twoways")}).
#' @param n.break Number of breaks.
#' If \code{n.break = 0}, it simply runs fixed effect model with shrinkage prior on coefficients.
#' @param mcmc MCMC iteration.
#' @param burn Burn-in period.
#' @param verbose Verbose.
#' @param thin Thinning.
#' @param c0 Hyperparam
#' @param d0 = 0.1
#' @param nu.shape =2.0
#' @param nu.rate =2.0
#' @param alpha  = 1
#'
#'
#' @author Jong Hee Park, and Soichiro Yamauchi \email{syamauchi@princeton.edu}
#'
#' @export
#'
#'
## Adaptive for each regime
adaptive.lasso <- function(y, x){
    fit.ridge <- cv.glmnet(y = y, x = x, type.measure="mse",
                           alpha=0, standardize = TRUE, family="gaussian")
    w3 <- 1/abs(matrix(coef(fit.ridge, s=fit.ridge$lambda.min)[, 1][2:(ncol(x)+1)] ))^1 ## Using gamma = 1
    w3[w3[,1] == Inf] <- 999999999 ## Replacing values estimated as Infinite for 999999999
    cv.adaptive <- cv.glmnet(x=x, y=y, family='gaussian', alpha=1, penalty.factor=w3)
    beta.adaptive <- coef(cv.adaptive)[-1]
    return(beta.adaptive)
}
