#' Detect a break number using different metrics 
#'
#'
#' @param Y Reponse tensor 
#' @param R Dimension of latent space. The default is 2. 
#' @param break.upper Upper threshold for break number detection.
#'  The default is \code{break.upper = 3}.
#'
#' @param burnin The number of burn-in iterations for the sampler.
#'
#' @param mcmc The number of MCMC iterations after burnin.
#'
#' @param thin The thinning interval used in the simulation.  The number of
#' MCMC iterations must be divisible by this value.
#'
#' @param verbose A switch which determines whether or not the progress of the
#' sampler is printed to the screen.  If \code{verbose} is greater than 0 the
#' iteration number, the \eqn{\beta} vector, and the error variance are
#' printed to the screen every \code{verbose}th iteration.
#'
#' 
#' @param degree.normal	A null model for degree correction. Users can choose "NULL", "eigen" or "Lsym."
#' "NULL" is no degree correction. "eigen" is a principal eigen-matrix consisting of
#' the first eigenvalue and the corresponding eigenvector. "
#' Lsym" is a modularity matrix. Default is "eigen."
#'
#' @param UL.Normal Transformation of sampled U. Users can choose "NULL", "Normal" or "Orthonormal."
#' "NULL" is no normalization. "Normal" is the standard normalization.
#' "Orthonormal" is the Gram-Schmidt orthgonalization. Default is "NULL."
#'
#' 
#' @param v0 \eqn{v_0/2} is the shape parameter for the inverse
#' Gamma prior on variance parameters for V.
#' If \code{v0 = NULL}, a value is computed from a test run of \code{NetworkStatic}.
#' 
#' @param v1 \eqn{v_1/2} is the scale parameter for the
#' inverse Gamma prior on variance parameters for V.
#' If \code{v1 = NULL}, a value is computed from a test run of \code{NetworkStatic}.
#'
#'
#' @export
#'
#'
#' @examples
#'    \dontrun{
#'    set.seed(1973)
#'    ## One break test
#'    out <- BridgeChangeSim(ntime=20, predictor = 10, n.break=1, constant.p =0, varying.p = 0.4, dgp.only=TRUE)
#'
#'    ## Fit multiple models for break number detection using Bayesian model comparison
#'    detect <- BreakDiagnostic(y=out$y.c, X=out$x.c)
#'    
#'    ## Look at the graph
#'    detect[[1]]; print(detect[[2]])
#'
#'    ## Two break test
#'    out <- BridgeChangeSim(ntime=20, predictor = 10, n.break=2, constant.p =0, varying.p = 0.4, dgp.only=TRUE)
#'
#'    ## Fit multiple models for break number detection using Bayesian model comparison
#'    detect <- BreakDiagnostic(y=out$y.c, X=out$x.c)
#'    
#'    ## Look at the graph
#'    detect[[1]]; print(detect[[2]])
#'   
#' }
#'
#' 

BreakDiagnostic <- function(y, X, mcmc=100, burn=100, verbose=100, thin=1, break.upper = 3){
    ## set.seed(11173)
    ## prior estimate
     
    ## model fit
    out <- as.list(rep(NA, break.upper))
    for(m in 1:(break.upper+1)){
        ## to save time and to be more conservative, use randomly generated initial states
        out[[m]] <- BridgeChangeReg(y = y, X = X, n.break = m-1,
                                    mcmc=mcmc, burn=burn, thin=thin, verbose=verbose,
                                    Waic = TRUE, marginal = TRUE)
    }
    
    ## diagnostic info
    Waic.holder <- marginal.holder <- loglike.holder <- rep(NA, 4)
    for(i in 1:(break.upper+1)){
        loglike.holder[i] <- attr(out[[i]], "loglike")
        marginal.holder[i] <- attr(out[[i]], "logmarglike")
        Waic.holder[i] <- attr(out[[i]], "Waic.out")[1]
    }
    ## loss
    loss.input <- out[-1]
    loss.out <- BreakPointLoss(loss.input, display=FALSE)[[1]]

    ## save model diagnostics
    result <- list("LogMarginal" = marginal.holder,
                   "Loglike" = loglike.holder,
                   "WAIC" = Waic.holder,
                   "Average Loss" = loss.out)

    test.curve1 <- -2*matrix(result[[1]], 1, break.upper +1, byrow=TRUE)
    test.curve2 <- -2*matrix(result[[2]], 1, break.upper +1, byrow=TRUE)
    test.curve3 <- matrix(result[[3]], 1, break.upper +1, byrow=TRUE)
    test.curve4 <- matrix(c(NA, result[[4]]), 1, break.upper +1, byrow=TRUE)

    test.curve <- rbind(test.curve1, test.curve2, test.curve3, test.curve4)
    test.curve <- data.frame(test.curve)
    colnames(test.curve) <- paste0("break", 0:break.upper)
    test.curve$Metric <- c("-2*LogMarginal", "-2*Loglike", "WAIC","Average Loss")
    data_long <- gather(test.curve, model, value, 
                        paste0("break", 0:break.upper), factor_key=TRUE)

    g1 <- ggplot(data= data_long, mapping = aes(x = model, y = value, group = Metric, color = Metric)) +
        geom_line(size=0.2) + geom_point(cex=3, alpha=1/2) + facet_wrap(~Metric, nrow=1, ncol=4, scales = "free_y") + 
        labs(x = "Model", y = "Value") + theme_bw() +
        theme(legend.position="none",
              legend.key = element_blank(),
              plot.title = element_text(hjust = 0.5))
    
    return(list(graph=g1, result=result))
}
