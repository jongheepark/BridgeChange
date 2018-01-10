######################################################################################################
## The idea is to develop a sparsity-induced prior-posterior model that fits data
## Bridge regression using mixture of normals representation.
## by JHP "Wed Oct 19 15:13:47 2016"
######################################################################################################

#' Sparse Change Point Model with Fixed Effect
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
BridgeFixedPanel <- function(formula, data, index, model, effect,
                             n.break = 1, mcmc = 100, burn = 100, verbose = 100, thin = 1,
                             b0=0, B0=1, c0 = 0.1, d0 = 0.1, r0 =  1, R0 = 1,
                             nu.shape = 2.0, nu.rate = 2.0, alpha = 1,
                             Waic = FALSE, marginal = FALSE) {
    call <- match.call()
    a = NULL; b = NULL
    ## ---------------------------------------------------- ##
    ## use plm package here
    ## transform data int pdata.frame object
    ## ---------------------------------------------------- ##

    pdata    <- pdata.frame(data, index)
    pformula <- pFormula(formula)


    ## get transformed Y and X
    X <- plm:::model.matrix.pFormula(formula, pdata, rhs = 1, model = model, effect = effect)
    y <- plm:::pmodel.response.pFormula(formula, pdata, model = model, effect = effect)


    ##
    ## centering X and Y?
    ##
    ## m <- n.break
    X <- apply(X, 2, scale)
    y <- scale(y)
    W <- matrix(0, length(y), 1)
    # if(inter){
    #     x1.1 <- data.frame(X)
    #     var.names <- colnames(X)
    #     x1.2 <- matrix(t(apply(x1.1, 1, combn, 2, prod)), nrow = nrow(X))
    #     newX <- as.matrix(cbind(x1.1, x1.2))
    #     colnames(newX) <- c(var.names, combn(var.names, 2, paste, collapse="-"))
    #     X <- newX
    # }
    data2 <- cbind(data[, index], y, X)
    var.names <- colnames(data2)
    colnames(data2) <- gsub("X", "", var.names)
    pdata2    <- pdata.frame(data2, index)
    pm <- plm(y~X-1, data = pdata2, model = model, effect = effect)

    ## plot
    coplot(pm$residuals ~ pdata[,index[2]]|pdata[,index[1]], data=pdata, number=num,
           overlap=.1, col="brown", type="l", panel = panel.smooth, xlab="panel residuals over group and time")
    ## coplot(pm$residuals ~ pdata[,index[2]]|pdata[,index[1]], type="b", data=pdata)


    subject.id <- as.numeric(as.factor(data[,index[1]]))
    time.id    <- as.numeric(as.factor(data[,index[2]]))

    ## ---------------------------------------------------- ##
    ## run change point model on demeaned data
    ## ---------------------------------------------------- ##
    output <- BridgeMixedPanel(subject.id = subject.id, time.id = time.id, y=y, X=X, W=W,
                               n.break = n.break, b0=b0, B0=B0, c0=c0, d0=d0, r0=r0, R0=R0,
                               ## inter should be always turned off because we already check it in the above!!
                               mcmc = mcmc, burn = burn, thin = thin, verbose=verbose, inter=FALSE,
                               nu.shape = 2.0, nu.rate = 2.0, alpha = 1, Waic = Waic, marginal = marginal, fixed = TRUE)

    attr(output, "title") <- "SparseChangeFixedPanel Posterior Sample"
    attr(output, "m")       <- n.break
    attr(output, "plm")     <- pm
    return(output)
}
