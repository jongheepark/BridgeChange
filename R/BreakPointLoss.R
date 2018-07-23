## ---------------------------------------------------- ##
## ---------------------------------------------------- ##
## The idea is to compute the average loss of sampled states from the expected break points
## by JHP "Fri Jul 20 09:57:45 2018"
## ---------------------------------------------------- ##
## ---------------------------------------------------- ##

mse <- function(x, mu){mean((x - mu)^2)}

trace.break <- function(tau.samp, tau){
  m <- length(tau)
  ## par(mfrow=c(1, m))
  if(m>1){
    ## for(i in 1:m){
    ##   plot(tau.samp[i,], type="l"); abline(h=tau[i], col="red")
    ## }
    out <- sapply(1:m, function(i){mse(tau.samp[i,], tau[i])})
  }else{
    ## plot(tau.samp, type="l"); abline(h=tau, col="red")
    out <- mse(tau.samp, tau)
  }
  return(out)
}


findBreakPoint <- function (mcmcout, start = 1) 
{
    out <- attr(mcmcout, "prob.state")
    
    y <- attr(mcmcout, "y")
    m <- attr(mcmcout, "m")
    ## time <- length(y)
    ## out <- sapply(1:(m+1), function(j){sapply(1:time, function(t){mean(sout[,t]==j)})})
    
    if (!is.ts(y)) 
        y <- ts(y, start)
    time.frame <- as.vector(time(y))
    if (m == 1) {
        pr.st <- c(0, diff(out[, (m + 1)]))
        pr.st[pr.st < 0] <- 0
        cp <- which(cumsum(pr.st) > 0.5)[1] - 1
    }else {
        cp <- rep(NA, m)
        for (i in 2:m) {
            pr.st <- c(0, diff(out[, i]))
            pr.st <- ifelse(pr.st < 0, 0, pr.st)
            cp[i - 1] <- which(cumsum(pr.st) > 0.5)[1] - 1
        }
        pr.st <- c(0, diff(out[, (m + 1)]))
        pr.st[pr.st < 0] <- 0
        cp[m] <- which(cumsum(pr.st) > 0.5)[1] - 1
    }
    if(sum(is.na(cp))>0){
        cat("\n At break = ", m, " one state is dominated by other states and a break point is not defined for this state. \n")
    }
    ## cp.means <- rep(NA, m + 1)
    ## cp.start <- c(1, cp + 1)
    ## cp.end <- c(cp, length(y))
    return(cp + time.frame[1])
}


#' Compute the Average Loss of Hidden State Changes from Expected Break Points
#'
#'
#' @param ... MCMC output objects. These have to be of class
#'   \code{mcmc} and have a \code{logmarglike} attribute. In what
#'   follows, we let \code{M} denote the total number of models to be
#'   compared.
#'
#' @param marginal If \code{marginal} is TRUE, \code{logmarglike} will be reported.
#'
#' @param display If \code{display} is TRUE, a plot of \code{ave.loss} will be produced. 
#' 
#' \code{BreakPointLoss}. ave.loss, logmarglike, State, Tau, Tau.samp
#' @return \code{BreakPointLoss} returns five objects. They are: \code{ave.loss} the expected loss for each model
#'   computed by the mean sqaured distance of hidden state changes from the expected break points
#' \textrm{Average Loss} = \frac{1}{M}\sum_{m=1}^{M}\left(\frac{1}{G}\sum_{g=1}^{G} (\bar{\tau}_m - \tau_{m}^{(g)})^2 \right);
#'   \code{logmarglike} the natural log of the marginal likelihood for each model; \code{State} sampled state vectors;
#'   \code{Tau} expected break points for each model; and \code{Tau.samp} sampled break points from hidden state draws.
#'
#' @export
#'
#'
#' @examples
#' \dontrun{
#' set.seed(1119)
#' n <- 99
#' ns <- 3
#' x1 <- runif(n)
#' true.beta1 <- c(2, -2)
#' true.beta2 <- c(-2,  -2)
#' true.beta3 <- c(0,  2)
#' true.Sigma <- c(1, 2)
#' true.s <- rep(1:ns, each=n/ns)
#' 
#' mu1 <- cbind(1, x1[true.s==1])%*%true.beta1
#' mu2 <- cbind(1, x1[true.s==2])%*%true.beta2
#' mu3 <- cbind(1, x1[true.s==3])%*%true.beta3
#' 
#' y <- as.ts(c(rnorm(n/ns, mu1, sd=sqrt(true.Sigma[1])), 
#'             rnorm(n/ns, mu2, sd=sqrt(true.Sigma[2])), 
#'              rnorm(n/ns, mu3, sd=sqrt(true.Sigma[1]))))
#' formula=y ~ x1
#' ## prior
#' b0 <- 0; B0 <- 1
#' sigma.mu=sd(y)
#' sigma.var=var(y)
#' mcmc = 1000
#' ## models
#' model1 <-  MCMCregressChange(formula, m=1, b0=b0, B0=B0, mcmc=mcmc, burnin=mcmc,
#'           sigma.mu=sigma.mu, sigma.var=sigma.var, marginal.likelihood="Chib95")
#' model2 <-  MCMCregressChange(formula, m=2, b0=b0, B0=B0, mcmc=mcmc, burnin=mcmc,
#'            sigma.mu=sigma.mu, sigma.var=sigma.var, marginal.likelihood="Chib95")
#' model3 <-  MCMCregressChange(formula, m=3, b0=b0, B0=B0, mcmc=mcmc, burnin=mcmc,
#'            sigma.mu=sigma.mu, sigma.var=sigma.var, marginal.likelihood="Chib95")
#' model4 <-  MCMCregressChange(formula, m=4, b0=b0, B0=B0, mcmc=mcmc, burnin=mcmc,
#'            sigma.mu=sigma.mu, sigma.var=sigma.var, marginal.likelihood="Chib95")
#' model5 <-  MCMCregressChange(formula, m=5, b0=b0, B0=B0, mcmc=mcmc, burnin=mcmc,
#'            sigma.mu=sigma.mu, sigma.var=sigma.var, marginal.likelihood="Chib95")
#'
#' out <- BreakPointLoss(model1, model2, model3, model4, model5, marginal=TRUE)
#'
#' print(out[["ave.loss"]])
#' }
#'
#' 
BreakPointLoss <- function(model.list, marginal=FALSE, display=TRUE){

    M <- length(model.list)
    this.call <- match.call()
    this.call.string <- deparse(this.call)
    this.call.string <- strsplit(this.call.string, "BreakPointLoss\\(")
    this.call.string <- this.call.string[[1]][length(this.call.string[[1]])]
    this.call.string <- strsplit(this.call.string, ",")
    
    break.number <- rep(NA, M)
    for (i in 1:M) {
        break.number[i] <- attr(model.list[[i]], "m")
        ## print(break.number[i])
        if(break.number[i] < 1){
            stop("no break model must be dropped\n")
        }
    }

    model.names <- paste0("break ", break.number)## c(model.names, this.call.string[[1]][i])
    ## If marginal, report marginal likelihood
    logmarglike <- NULL
    if(marginal){
        logmarglike <- BayesFactor(...)[[3]]
    }
    State <- Tau <- Tau.samp <- as.list(rep(NA, M))

    for (i in 1:M) {
        State[[i]] <- attr(model.list[[i]], "s.store")
        Tau[[i]] <- findBreakPoint(model.list[[i]])
        Tau.samp[[i]] <- sapply(1:nrow(State[[i]]), function(j){which(diff(State[[i]][j,])==1)+1})
    }
    
    ## Report Average Loss
    ave.loss <- rep(NA, M)
    for (i in 1:M) {
        ave.loss[i] <- mean(trace.break(Tau.samp[[i]], Tau[[i]]))
    }

    if(display){
        plot(ave.loss, type="o", xlab="Model", ylab="Loss", main="Average Loss", 
             axes=FALSE)
        axis(1, at=break.number, labels=model.names); axis(2)
        abline(v = which.min(ave.loss), lty=3, col="red")
    }
    return(list(ave.loss = ave.loss, logmarglike=logmarglike, State=State, Tau=Tau, Tau.samp=Tau.samp))
}
