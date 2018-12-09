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


BridgeFixedPanelHybrid <- function(formula, data, index, model, effect,
                             standardize = TRUE, inter = FALSE, 
                             n.break = 1, 
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
    
    if(inter){
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
        dat.sd <- c(ysd, Xsd)
        X <- scale(X)
        y <- scale(as.vector(y))        
    }


    subject.id <- as.numeric(as.factor(data[,index[1]]))
    time.id    <- as.numeric(as.factor(data[,index[2]]))

    ## ---------------------------------------------------- ##
    ## run change point model on demeaned data
    ## ---------------------------------------------------- ##
    output <- BridgeMixedPanel(subject.id = subject.id, time.id = time.id, y=as.vector(y), X=X, W=W,
                               n.break = n.break, b0=b0, B0=B0, c0=c0, d0=d0, r0=r0, R0=R0,
                               standardize = FALSE,
                               mcmc = mcmc, burn = burn, thin = thin, verbose=verbose, 
                               nu.shape = 2.0, nu.rate = 2.0, alpha = 1, Waic = Waic, marginal = marginal, fixed = TRUE,
                               unscaled.Y = unscaled.Y, unscaled.X = unscaled.X)

    ## data preparation for regime-wise regression
    state <- round(apply(attr(output, "s.store"), 2, mean))
    unique.time.index <- sort(unique(attr(plmX,"index")[,2]))
    n.state <- length(unique(state))
    y.list <- x.list <- as.list(rep(NA, n.state))
    for(i in 1:n.state){
        x.list[[i]] <- X[ is.element(attr(plmX,"index")[,2], unique.time.index[state==i]), ]
        y.list[[i]] <- y[ is.element(attr(plmX,"index")[,2], unique.time.index[state==i])]    
    }
    ## Variable selection using adaptive lasso
    res.inter <- sapply(1:n.state, function(i){adaptive.lasso(y.list[[i]], x.list[[i]])})
    rownames(res.inter) <- colnames(X)
    colnames(res.inter) <- paste0("Regime", 1:n.state)
    
    attr(output, "hybrid") <- res.inter
    attr(output, "title")  <- "SparseChangeFixedPanel Posterior Sample"
    attr(output, "m")      <- n.break
    if(standardize) attr(output, "dat.sd") <- dat.sd
    # attr(output, "plm")     <- pm
    # class(output) <- c("mcmc", "BridgeChange")
    return(output)
}
