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
#' @param interaction If interaction = 1, no interaciton. If interaction = 2, only two-way interaciton. Interaction can be up to K, which is the rank of the model matrix. 
#' @param n.break Number of breaks.
#' If \code{n.break = 0}, it simply runs fixed effect model with shrinkage prior on coefficients.
#' @param burnin The number of burn-in iterations for the sampler.
#'
#' @param mcmc The number of MCMC iterations after burn-in.
#'
#' @param thin The thinning interval used in the simulation.  The
#'   number of MCMC iterations must be divisible by this value.
#'
#' @param verbose A switch which determines whether or not the
#'   progress of the sampler is printed to the screen.  If
#'   \code{verbose} is greater than 0, the iteration number and the
#'   posterior density samples are printed to the screen every
#'   \code{verbose}th iteration.
#'
#' @param alpha.MH  If \code{TRUE}, alpha is updated by MH algorithm. By default, it is updated by griddy gibbs.
#' @param c0 Hyperparam
#' @param d0 = 0.1
#' @param nu.shape =2.0
#' @param nu.rate =2.0
#' @param alpha  = 1
#' @param a \eqn{a} is the shape1 beta prior for transition
#'   probabilities.  By default, the expected duration is computed and
#'   corresponding a and b values are assigned. The expected duration
#'   is the sample period divided by the number of states.
#'
#' @param b \eqn{b} is the shape2 beta prior for transition
#'   probabilities.  By default, the expected duration is computed and
#'   corresponding a and b values are assigned. The expected duration
#'   is the sample period divided by the number of states.

#'
#'
#' @author Jong Hee Park, and Soichiro Yamauchi \email{syamauchi@princeton.edu}
#' @importFrom plm pmodel.response
#' @importFrom utils getFromNamespace
#' @example examples/fixed_panel_eg.R

#' 
#' @export
#' 
BridgeFixedPanel <- function(formula, data, index, model, effect,
                             standardize = TRUE,
                             interaction = 1,
                             a = NULL, b = NULL,
                             n.break = 1, alpha.MH = FALSE,
                             mcmc = 100, burn = 100, verbose = 100, thin = 1,
                             c0 = 0.1, d0 = 0.1, r0 =  1, R0 = 1,
                             nu.shape = 2.0, nu.rate = 2.0, alpha = 1,
                             Waic = FALSE, marginal = FALSE) {
    call <- match.call()
    
    
    ## @importMethodsFrom does not work 
    model.matrix <- getFromNamespace("model.matrix.pFormula", "plm")
    
    ## ---------------------------------------------------- ##
    ## use plm package here
    ## transform data int pdata.frame object
    ## ---------------------------------------------------- ##
    # if (standardize) {
    #     dat.sd <- apply(data[, !(colnames(data) %in% index)], 2, sd)
    #     data <- data.frame(cbind(scale(data[, !(colnames(data) %in% index)]), data[,index]))
    # }
    
    suppressWarnings(pdata    <- pdata.frame(data, index))
    suppressWarnings(pformula <- pFormula(formula))
    
    suppressWarnings(X <- plm:::model.matrix.pFormula(pformula, pdata, rhs = 1, model = model, effect = effect))
    suppressWarnings(y <- plm:::pmodel.response(pformula, pdata, model = model, effect = effect))

    ## Drop covariates with all zero
    if(sum(apply(X, 2, sd) == 0)>0){
        cat("Some interactions (", sum(apply(X, 2, sd) == 0) , ") are all zero. So they are removed!\n")
        X <- X[, apply(X, 2, sd)!=0]
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
    ##
    ## centering X and Y?
    ##
    m <- n.break
    #     
    
    W <- matrix(0, length(y), 1)
    
    interaction <- min(ncol(X), interaction)
    if(interaction>1){
        newX <- list()
        newX[[1]] <- X
        x1.1 <- data.frame(X)
        var.names <- colnames(X)
        for(j in 2:interaction){
            x1.2 <- matrix(t(apply(x1.1, 1, combn, j, prod)), nrow = nrow(X))
            newX[[j]] <- as.matrix(x1.2)
            colnames(newX[[j]]) <- c(combn(var.names, j, paste, collapse="-"))
        }
        X <- Reduce(cbind, newX)
        
        ## Drop covariates with all zero
        if(sum(apply(X, 2, sd) == 0)>0){
            cat("Some interactions (", sum(apply(X, 2, sd) == 0) , ") are all zero. So they are removed!\n")
            X <- X[, apply(X, 2, sd)!=0]
        }
    }

    subject.id <- as.numeric(as.factor(data[,index[1]]))
    time.id    <- as.numeric(as.factor(data[,index[2]]))

    ## ---------------------------------------------------- ##
    ## run change point model on demeaned data
    ## ---------------------------------------------------- ##
    output <- BridgeMixedPanel(subject.id = subject.id, time.id = time.id, y=as.vector(y), X=X, W=W,
                               n.break = n.break, c0=c0, d0=d0, r0=r0, R0=R0,
                               standardize = FALSE, alpha.MH = alpha.MH, 
                               mcmc = mcmc, burn = burn, thin = thin, verbose=verbose, 
                               nu.shape = 2.0, nu.rate = 2.0, alpha = 1, Waic = Waic,
                               marginal = marginal, fixed = TRUE,
                               unscaled.Y = unscaled.Y, unscaled.X = unscaled.X)

    attr(output, "title")  <- "BridgeChangeFixedPanel Posterior Sample"
    attr(output, "m")      <- n.break
    if(standardize) attr(output, "dat.sd") <- dat.sd
    # attr(output, "plm")     <- pm
    # class(output) <- c("mcmc", "BridgeChange")
    return(output)
}
