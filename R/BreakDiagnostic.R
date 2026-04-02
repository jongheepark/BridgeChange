#' Detect a break number using different metrics
#'
#' Fits models with 0 to \code{break.upper} breaks and compares them using
#' log marginal likelihood, log likelihood, WAIC, and average loss.
#'
#' @param y Response vector.
#' @param X Design matrix.
#' @param mcmc The number of MCMC iterations after burn-in. Default is 100.
#' @param burn The number of burn-in iterations. Default is 100.
#' @param verbose Print progress every \code{verbose}th iteration. Default is 100.
#' @param thin Thinning interval. Default is 1.
#' @param break.upper Upper limit for break number detection. Default is 3.
#'
#' @return A list containing a ggplot diagnostic graph and a list of model comparison metrics.
#'
#' @export
#'
#'
#' @examples
#'    \dontrun{
#'    set.seed(1973)
#'    ## One break test
#'    out <- BridgeChangeSim(ntime=20, predictor = 10,
#'      n.break=1, constant.p=0, varying.p=0.4, dgp.only=TRUE)
#'
#'    ## Fit multiple models for break number detection
#'    detect <- BreakDiagnostic(y=out$y.c, X=out$x.c)
#'
#'    ## Look at the graph
#'    detect[[1]]; print(detect[[2]])
#'
#'    ## Two break test
#'    out <- BridgeChangeSim(ntime=20, predictor = 10,
#'      n.break=2, constant.p=0, varying.p=0.4, dgp.only=TRUE)
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
        out[[m]] <- BridgeChangeRegHybrid(y = y, X = X, n.break = m-1,
                                    mcmc=mcmc, burn=burn, thin=thin, verbose=verbose,
                                    waic = TRUE, marginal = TRUE)
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
    break_cols <- paste0("break", 0:break.upper)
    colnames(test.curve) <- break_cols
    test.curve$Metric <- c("-2*LogMarginal", "-2*Loglike", "WAIC","Average Loss")

    ## reshape to long format (replaces tidyr::gather)
    data_long <- data.frame(
        Metric = rep(test.curve$Metric, each = length(break_cols)),
        model = factor(rep(break_cols, times = nrow(test.curve)), levels = break_cols),
        value = as.vector(t(as.matrix(test.curve[, break_cols])))
    )

    g1 <- ggplot(data= data_long, mapping = aes(x = model, y = value, group = Metric, color = Metric)) +
        geom_line(linewidth=0.2) + geom_point(size=3, alpha=1/2) + facet_wrap(~Metric, nrow=1, ncol=4, scales = "free_y") +
        labs(x = "Model", y = "Value") + theme_bw() +
        theme(legend.position="none",
              legend.key = element_blank(),
              plot.title = element_text(hjust = 0.5))
    
    return(list(graph=g1, result=result))
}
