## From BayesLogit
################################################################################
## PG(1.0, Z) - FOLLOWING DEVROYE ##
################################################################################

TRUNC = 0.64
cutoff = 1 / TRUNC;
## Sample from PG(n, Z) using Devroye-like method.
## n is a natural number and z is a positive real.
##------------------------------------------------------------------------------
rpg.devroye <- function(num=1, n=1, z=0.0)
{
  z = array(z, num);
  n = array(n, num);

  total.trials = 0;

  x = rep(0, num);
  for (i in 1:num) {
    x[i] = 0;
    for (j in 1:n[i]) {
      ## x[i] = x[i] + rpg.devroye.1(z[i])
      temp = rpg.devroye.1(z[i]);
      x[i] = x[i] + temp$x;
      total.trials = total.trials + temp$n;
    }
  }
  x
  list("x"=x, "rate"=sum(n)/total.trials)
}
## Samples from PG(n=1.0, psi=Z)
## ------------------------------------------------------------------------------
rpg.devroye.1 <- function(Z)
{
  Z = abs(Z) * 0.5;

  ## PG(1,z) = 1/4 J*(1,Z/2)
  fz = pi^2 / 8 + Z^2 / 2;
  ## p = (0.5 * pi) * exp( -1.0 * fz * TRUNC) / fz;
  ## q = 2 * exp(-1.0 * Z) * pigauss(TRUNC, 1.0/Z, 1.0);

  num.trials = 0;
  total.iter = 0;

  while (TRUE)
    {
      num.trials = num.trials + 1;

      if ( runif(1) < mass.texpon(Z) ) {
        ## Truncated Exponential
        X = TRUNC + rexp(1) / fz
      }
      else {
        ## Truncated Inverse Normal
        X = rtigauss(Z)
      }

      ## C = cosh(Z) * exp( -0.5 * Z^2 * X )

      ## Don't need to multiply everything by C, since it cancels in inequality.
      S = a.coef(0,X)
      Y = runif(1)*S
      n = 0

      while (TRUE)
        {
          n = n + 1
          total.iter = total.iter + 1;
          if ( n %% 2 == 1 )
            {
              S = S - a.coef(n,X)
              if ( Y<=S ) break
            }
          else
            {
              S = S + a.coef(n,X)
              if ( Y>S ) break
            }
        }

      if ( Y<=S ) break
    }

  ## 0.25 * X
  list("x"=0.25 * X, "n"=num.trials, "total.iter"=total.iter)
}
## rtigauss - sample from truncated Inv-Gauss(1/abs(Z), 1.0) 1_{(0, TRUNC)}.
##------------------------------------------------------------------------------
rtigauss <- function(Z, R=TRUNC)
{
  Z = abs(Z);
  mu = 1/Z;
  X = R + 1;
  if (mu > R) {
    alpha = 0.0;
    while (runif(1) > alpha) {
      ## X = R + 1
      ## while (X > R) {
      ##     X = 1.0 / rgamma(1, 0.5, rate=0.5);
      ## }
      E = rexp(2)
      while ( E[1]^2 > 2 * E[2] / R) {
        E = rexp(2)
      }
      X = R / (1 + R*E[1])^2
      alpha = exp(-0.5 * Z^2 * X);
    }
  }
  else {
    while (X > R) {
      lambda = 1.0;
      Y = rnorm(1)^2;
      X = mu + 0.5 * mu^2 / lambda * Y -
        0.5 * mu / lambda * sqrt(4 * mu * lambda * Y + (mu * Y)^2);
      if ( runif(1) > mu / (mu + X) ) {
        X = mu^2 / X;
      }
    }
  }
  X;
}
## Calculate coefficient n in density of PG(1.0, 0.0), i.e. J* from Devroye.
##------------------------------------------------------------------------------
a.coef <- function(n,x)
{
  if ( x>TRUNC )
    pi * (n+0.5) * exp( -(n+0.5)^2*pi^2*x/2 )
  else
    (2/pi/x)^1.5 * pi * (n+0.5) * exp( -2*(n+0.5)^2/x )
}

## when z hits the Infs, use 300 or -300 instead!
"z.latent.internal" <-function(Y, mu, j){
    muj <-  mu[j]
    pj <-  pnorm(-muj)
    yj  <-  Y[j]
    ## sdj <-  sd[j]
    uj  <-  runif(1)
    z   <-  ifelse(yj==0, muj+qnorm(uj*pj), muj+qnorm(pj+uj*(1-pj)))
    ## If pj is 1, qnorm(pj+uj*(1-pj)) becomes Inf
    ## So, use muj instead of new z if pj is 1.
    if (z==-Inf){
        out <- -300
    } else if (z==Inf){
        out <- 300
    } else {
      out <- z
    }
  return(out)
}
mass.texpon <- function(Z)
{
  x = TRUNC;
  fz = pi^2 / 8 + Z^2 / 2;
  b = sqrt(1.0 / x) * (x * Z - 1);
  a = -1.0 * sqrt(1.0 / x) * (x * Z + 1);

  x0 = log(fz) + fz * TRUNC;
  xb = x0 - Z + pnorm(b, log.p=TRUE);
  xa = x0 + Z + pnorm(a, log.p=TRUE);

  qdivp = 4 / pi * ( exp(xb) + exp(xa) );

  1.0 / (1.0 + qdivp);
}


"z.latent.internal2" <-function(Y, mu){
  N <- length(Y)
  left  = rep(0, N);
  left [Y==0]=-Inf;
  right = rep(0, N);
  right[Y==1]= Inf;
  z = msm::rtnorm(N, mu, 1.0, left, right);
  return(z)
}

## function to generate all interaction terms
interFun <- function(formula, data){
    lm1 <- lm(formula, data=data, x=TRUE, y=TRUE)
    X <- lm1$x[, -1]
    if(is.null(colnames(data))){
        var.names <- paste0("x", 1:ncol(X))
    }else{
        var.names <- colnames(X)
    }
    x1.1 <- data.frame(X)
    x1.2 <- t(apply(x1.1, 1, combn, 2, prod))
    x1.3 <- as.matrix(cbind(x1.1, x1.2))
    colnames(x1.3) <- c(var.names, combn(var.names, 2, paste, collapse="-"))
    out <- cbind(lm1$y, x1.3)
    colnames(out)[1] <- "y"
    return(out)
}
## centering required
centerdata <- function(X, all=TRUE){
    new.X <- X
    K <- dim(X)[2]
    col.mean <- colMeans(X)
    col.sd <- apply(X, 2, sd)
    for(k in 1:K){
        if(!all){
            ## if binary data should be untouched
            if(length(unique(X[,k])) == 2){ new.X[,k] <- X[,k] }
        }
        new.X[,k] <- (X[,k] - col.mean[k]) ## /(col.sd[k])
    }
    return(new.X)
}
center <- function(x){
    out <- (x - mean(x, na.rm=TRUE))/sd(x, na.rm=TRUE)
}
within.transform <- function(y, x, generation){
    newx <- x
    newy <- y
    for(i in unique(generation)){
        sub <- x[generation == i,]
        suby <- y[generation == i]
        dx <- colMeans(sub)
        newx[generation == i,] <- sweep(sub, 2, dx, FUN="-")
        newy[generation == i] <- suby - mean(suby, na.rm=TRUE)
    }
    return(list(newy, newx))
}
MSEcompare <- function(outlist){
    N.model <- length(outlist)
    breaks <- lapply(outlist, attr, "m")
    outl <- lapply(outlist, attr, "MSE")
    out <- matrix(outl, 1, N.model)
    colnames(out) <- paste0("break ", breaks)
    return(out)
}

#' Compare Models based on WAIC 
#' @param output list of \code{\link{BridgeChange}} objects, whch are typically outputs from \code{BridgeChangeReg}.
#' @param print boolean; if \code{TRUE} selected model is printed.
#' @return A vector of WAIC
WaicCompare <- function(outlist, print = TRUE){
    N.model <- length(outlist)
    breaks <- lapply(outlist, attr, "m")
    outl <- lapply(outlist, attr, "Waic.out")
    outm <- matrix(unlist(outl), N.model, 8, byrow=TRUE)
    out <- matrix(outm[,1], 1, N.model)
    colnames(out) <- paste0("break ", breaks)
    
    if(isTRUE(print)) {
      cat("\n")
      cat("Selected model = break", unlist(breaks)[which.min(out[1,])],'\n')        
      cat("\n")
      print.default(out[1,], quote = FALSE, digit = 5)
    }
    return(out)
}
MarginalCompare <- function(outlist){
    N.model <- length(outlist)
    breaks <- lapply(outlist, attr, "m")
    outl <- lapply(outlist, attr, "logmarglike")
    out <- matrix(outl, 1, N.model)
    colnames(out) <- paste0("break ", breaks)
    return(out)
}

## code by Gelman and Vehtari (2014)
colVars <- function(a) {
  n <- dim(a)[[1]]; c <- dim(a)[[2]];
  out <- .colMeans(((a - matrix(.colMeans(a, n, c),
    nrow = n, ncol = c, byrow = TRUE)) ^ 2), n, c) * n / (n - 1)
  return(out)
}

waic_calc <- function(log_lik){
  ## log_lik <- extract (stanfit, "log_lik")$log_lik
  dim(log_lik) <- if (length(dim(log_lik))==1) c(length(log_lik),1) else
  c(dim(log_lik)[1], prod(dim(log_lik)[2:length(dim(log_lik))]))
  S <- nrow(log_lik)
  n <- ncol(log_lik)
  lpd <- log(colMeans(exp(log_lik)))
  p_waic <- colVars(log_lik)
  p_waic1 <- 2*(log(colMeans(exp(log_lik))) - colMeans(log_lik))
  elpd_waic <- lpd - p_waic
  waic <- -2*elpd_waic
  waic2 <- -2*(lpd - p_waic1)
  loo_weights_raw <- 1/exp(log_lik-max(log_lik))
  loo_weights_normalized <- loo_weights_raw/
  matrix(colMeans(loo_weights_raw),nrow=S,ncol=n,byrow=TRUE)
  loo_weights_regularized <- pmin (loo_weights_normalized, sqrt(S))
  elpd_loo <- log(colMeans(exp(log_lik)*loo_weights_regularized)/
  colMeans(loo_weights_regularized))
  p_loo <- lpd - elpd_loo
  pointwise <- cbind(waic,waic2, lpd,p_waic,p_waic1, elpd_waic,p_loo,elpd_loo)
  total <- colSums(pointwise)
  ## this is strange? SE = s/sqrt(n)
  se <- sqrt(n*colVars(pointwise))
  ## se <- sqrt(colVars(pointwise)/n)
  waic.ci <- c(total[1] - 1.96*se[1], total[1] + 1.96*se[1])
  return(list(waic=total["waic"], waic2=total["waic2"],  elpd_waic=total["elpd_waic"],
  p_waic=total["p_waic"], p_waic1 = total["p_waic1"], elpd_loo=total["elpd_loo"], p_loo=total["p_loo"],
  pointwise=pointwise, total=total, se=se, waic.ci=waic.ci))
}
#

addTrans <- function(color,trans){
    ## thanks to Sacha Epskamp
    ## originally posted in
    ## http://stackoverflow.com/questions/12995683/any-way-to-make-plot-points-in-scatterplot-more-transparent-in-r
    ## This function adds transparancy to a color.
    ## Define transparancy with an integer between 0 and 255
    ## 0 being fully transparant and 255 being fully visable
    ## Works with either color and trans a vector of equal length,
    ## or one of the two of length 1.

    if (length(color)!=length(trans)&!any(c(length(color),length(trans))==1)) stop("Vector lengths not correct")
    if (length(color)==1 & length(trans)>1) color <- rep(color,length(trans))
    if (length(trans)==1 & length(color)>1) trans <- rep(trans,length(color))

    num2hex <- function(x)
        {
            hex <- unlist(strsplit("0123456789ABCDEF",split=""))
            return(paste(hex[(x-x%%16)/16+1],hex[x%%16+1],sep=""))
        }
    rgb <- rbind(col2rgb(color),trans)
    res <- paste("#",apply(apply(rgb,2,num2hex),2,paste,collapse=""),sep="")
    return(res)
}

DP.update <- function(F=F, state=state, switch1=switch1, a=a, b=b, time.of.change=time.of.change, j = j){
    ps.ss <- rep(NA, 2)
    recent.state <- state[time.of.change[j]- 1]   ## s[i]
    if(!is.na(time.of.change[j]) & (time.of.change[j] + 1 < nrow(F))){
        next.state <- state[time.of.change[j] + 1]  ## s[i + 1]
    }
    else{
        next.state <- recent.state + 1
    }
    nii <- switch1[recent.state, recent.state]
    ni1i1 <- switch1[next.state, next.state]
    ## merging prob.
    ps.ss[1] <- (b)/(a + b + nii) *
        ((ni1i1 + a)/(ni1i1 + a + b)) * F[time.of.change[j], next.state]
    ## splitting prob.
    ps.ss[2] <- (a + nii)/(a + b + nii) *
        (b/(nii + a + b + 1)) * F[time.of.change[j], recent.state]
    ## normalize
    pstyn   <-  ps.ss/sum(ps.ss)
    ## merge to next.state with prob of pstyn[1]. OW, stay at recent.state
    if(runif(1) < pstyn[1]){
        ## After sampling, if merge into t + 1, move to next one.
        new.state   <-  next.state
        cat("j ", j , " is changed from ", recent.state, " into ", next.state, "\n")
    }
    else{
        new.state <- recent.state
    }
    return(new.state)
}


## Dirichlet Process Prior CP estimate with Fixed a and b
DP.state.sampler <- function(state=state, m=m, y=y, X=X, a=3, b=2, beta=beta, sig2=sig2, P = P){
    T   <-  length(y)
    F   <-  matrix(NA, T, m+1)     # storage for the Filtered probabilities
    pr1 <-  c(1,rep(0, m))         # initial probability Pr(s=k|Y0, lambda)
    for (t in 1:T){
        mu <- X[t,] %*% t(beta)
        py  <-  sapply(c(1:(m+1)), function(i){
            dnorm(y[t], mean= mu[i], sd = sqrt(sig2[i]))})
        if(t==1) {
            pstyt1 = pr1
            F[t,]   <-  py/sum(py)
        }
        else {
            pstyt1     <- F[t-1,]%*%P
            unnorm.pstyt    <- pstyt1*py
            pstyt   <-  unnorm.pstyt/sum(unnorm.pstyt) # Pr(st|Yt)
            F[t,]   <-  pstyt
        }
    }

    new.s  <- state ## holder for new state variables
    ## ps     <-  matrix(NA, T, m+1) ## holder for state probabilities

    ##
    switch1 <- BridgeChange:::switchg(state, m=m)
    aba <- a/(b + a)
    bba <- b/(b + a)
    ps.ss <- rep(NA, 2)
    time.of.change <- which(diff(state) == 1)
    n.time.of.change <- length(time.of.change)

    ## case 1: for t = 1
    if(state[1] != state[2]){
        ## merging prob. s[1] == s[2]
        ps.ss[1] <- bba * (switch1[2,2] + a)/(switch1[2,2] + b + a) * F[1,1] ## eq 17. on p.281
        ## splitting prob. s[1] != s[2]
        ps.ss[2] <- aba * bba * F[1,2]
        ## normalize
        pstyn   <-  ps.ss/sum(ps.ss)
        ## merge to s[2] with prob of pstyn[1]. OW, stay at s[1]
        new.s[1]   <-  ifelse (runif(1) < pstyn[1], state[2], state[1])
    }
    if(state[T-1] != state[T]){
        ## merge to s[T-1]
        nTT <- switch1[state[T-1],state[T-1]]
        ps.ss[1] <- (nTT + a)/(nTT + b + a) * F[T,state[T-1]] ## eq 17. on p.281
        ## unchanged
        ps.ss[2] <- (b/(nTT + b + a))* F[T, state[T]]
        ## normalize
        pstyn   <-  ps.ss/sum(ps.ss)
        ## ps[T,] <-  pstyn
        ## merge to s[2] with prob of pstyn[1]. OW, stay at s[1]
        new.s[T]   <-  ifelse (runif(1) < pstyn[1], state[T-1], state[T])
    }
    ## case 2: for t = 2,,,T-1
    for(j in 1:n.time.of.change){
        if(time.of.change[j] != 1 & time.of.change[j] !=T){
            new.s[time.of.change[j]] <- DP.update(F=F, state=state,
                                                  switch1=switch1, a=a, b=b,
                                                  time.of.change=time.of.change, j = j)
        }
    }
    new.ss <- new.s
    ##################################################
    ## NOTE by JHP on "Mon Oct 24 11:37:53 2016"
    ##################################################
    ## This part of the code is problematic.
    ## They say we need to resample when new disjunctures of state (e.g. 1, 1, 2, 2 -> 1, 1, 1, 2) occur.
    ## That is, the last 2 must be resampled to see whether this needs to be down into 1 or jump into 2.
    ## Then, why we stop there? Let's say we sample new state to be 1. Then, 1,1,1,1 leads us to sample
    ## the next state, which is 1,1,1,1,2 because we don't know whether the last 2 is still valid given
    ## the new states of 1,1,1,1 sequence.....In this way, the resampling can't be stopped....

    ## when we stay, we resample t+1
    if(sum(state[time.of.change] == new.s[time.of.change])>0){
         time.of.stay <- time.of.change[which(state[time.of.change] == new.s[time.of.change])] + 1
         cat("--------------------------- \n")
         cat("state is sampled differently! \n")
         cat("time.of.change: ", time.of.change,  "\n")
         cat("time.of.stay: ", time.of.stay,  "\n")
         cat("time.of.change[which(state[time.of.change] == new.s[time.of.change])] + 1: ",
          time.of.change[which(state[time.of.change] == new.s[time.of.change])] + 1,  "\n")
         cat("old state: ", state,  "\n")
         cat("new state: ", new.s,  "\n")
          n.time.of.stay <- length(time.of.stay)
         for(j in 1:n.time.of.stay){
             new.ss[time.of.stay[j]] <- DP.update(F=F, state=new.s, switch1=switch1,
                                                  a=a, b=b, time.of.change=time.of.stay, j = j)
         }
     } ## end of if(length(s[time.of.change] ==new.s[time.of.change])>0){
    ## recent.state <- s[time.of.change[j]] ## s[i]
    ## next.state <- s[time.of.change[j+1]] ## s[i + 1]
    ## nii <- switch1[recent.state, recent.state]
    ## ni1i1 <- switch1[next.state, next.state]
    ## merging prob.
    ## ps.ss[1] <- (b)/(a + b + nii) *
    ##     ((ni1i1 + a)/(ni1i1 + a + b)) * F[time.of.change[j], next.state]
    ## splitting prob.
    ## ps.ss[2] <- (a + nii)/(a + b + nii) *
    ##     (b/(nii + a + b + 1)) * F[time.of.change[j], recent.state]
    ## normalize
    ## pstyn   <-  ps.ss/sum(ps.ss)
    ## merge to next.state with prob of pstyn[1]. OW, stay at recent.state
    ## if(runif(1) < pstyn[1]){
    ##     ## After sampling, if merge into t + 1, move to next one.
    ##     new.s[time.of.change[j]]   <-  next.state
    ## }
    ## else{
    ##    ## if stays, sample time.of.change + 1.
    ##     while(time.of.change[j] + 1 == time.of.change[j] + 1)

    ## case 3: for t = T
    ## merging prob.
    ## renormalize state vector
    ## new.sss <- as.numeric(as.factor(new.ss))
    length.of.each.state <- diff(c(0, c(which(diff(new.ss) !=0) + 1), T))
    if(sum(length.of.each.state == 0)>0){
        ## if zero count, add 0. zero count means ...., 10, 11, 10. The final transition is to another state...
        length.of.each.state[which(length.of.each.state == 0)] <- 0
    }
    new.sss <- rep(1:length(length.of.each.state),  length.of.each.state)

    ## name and report outputs
    out <-  list(new.sss)
    names(out) <- c("state")
    return(out)
}

###############################################
rdirichlet.cp <- function(n, alpha){
    ## pick n random deviates from the Dirichlet function
    ## with shape parameters alpha
    ## alpha can contain zero
    ## returns a matrix

    ## probability matrix for storage
    col<-length(alpha)
    out<-matrix(NA, n, col)

    ## find zero elements in alpha
    ## assign zero probability for alpha[alpha==0]
    out[,which(alpha==0)]<-0

    ## compute nonzero part and store them at dir.out
    a<-alpha[alpha!=0]
    l<-length(a);
    x<-matrix(rgamma(l*n,a),ncol=l,byrow=TRUE);
    sm<-x%*%rep(1,l);
    dir.out<-x/as.vector(sm);

    ## combine zero and nonzero parts prob
    out[,which(alpha!=0)]<-dir.out
    return(out)
}


###############################################
switchg <-
###############################################
    ## compute frequency of elements in state vector
    function(s1, m){
    ## for one step ahead case only
    ## example: s1<-c(1,1,1,2,2,2,3,4,4,4,4,5,5)
    ## switchg(s1)
        ## smax <- max(s1)
        out <- matrix(0, m+1 , m + 1)
        us <- unique(s1)
        tab.s <- table(s1)
        ## all P(i,i+1) are 1
        for (i in 1:length(us)){
            out[us[i],us[i+1]] <- 1
        }
        ## all P(i,i+1) are 1
        for (i in 1:length(us)){
            out[us[i],us[i]] <- tab.s[i]
        }
        return(out)
    }


###############################################
trans.mat.prior<-
###############################################
    ## function that generates transition prior (matrix)
    function(m, n, a=NULL, b=NULL){
    ## m is the number of changepoints
    ## n is the length of sample period
    ## output is m+1 by m+1 transition matrix
    ## (If a and b are not specified, this function
    ## automatically generates prior transition matrix
    ## that imposes a change point in the middle but
    ## user can change it by manipulating a and b)

        if (!is.null(a)|!is.null(b)){
            a <- a
            b <- b
        }
        else {expected.duration <- round(n/(m+1))
            b <- 0.1
            a <- b*expected.duration
        }

        ## make a transition matrix
        trans <- diag(m+1)

        # put a as diagonal elements except the last row
        diag(trans)[1:m]<-rep(a, m)

        # put b in trans[i, i+1]
        for (i in 1:m){trans[i, i+1]<-b}
        return(trans)
    }

llh.alpha <- function(alpha, s)
{
  p = length(s);
  p * log(alpha) - p * lgamma(1/alpha) - sum(exp(alpha * s))
}

## this is for alpha = (0,1) if alpha>1, dbeta(a.new, pr.a, pr.b, log=TRUE) = -Inf and log.accept = NaN
draw.alpha <- function(alpha, beta, tau, pr.a=1.0, pr.b=1.0, ep=0.2)
{
                                        # if (sum(abs(beta) <= 1e-5) == 0 & abs(tau) > 1e-5) {
    s = log(abs(beta / tau));
                                        # } else {
                                        #   beta[abs(beta) <= 1e-5] <- 1e-4
                                        #   if (abs(tau) <= 1e-5) tau <- 1e-4
                                        #   s = log(abs(beta / tau));
                                        # }

  a.old = alpha;
  lb = 0.0
  ub = 2.0

  ## We might want to change these bounds if we have some belief about the value
  ## of alpha.
  l.new = max(c(lb, a.old-ep))
  r.new = min(c(ub, a.old+ep))
  d.new = r.new - l.new
  a.new = runif(1, l.new, r.new);

  d.old = min(c(ub, a.new+ep)) -  max(c(lb, a.new-ep))

  if (ub < 1) {
      log.accept = llh.alpha(a.new, s) - llh.alpha(a.old, s) +
          dbeta(a.new, pr.a, pr.b, log=TRUE) - dbeta(a.old, pr.a, pr.b, log=TRUE) +
              log(d.old) - log(d.new);
  } else {
      log.accept = llh.alpha(a.new, s) - llh.alpha(a.old, s) +
          # dunif(a.new, lb, ub, log=TRUE) - dunif(a.old, lb, ub, log=TRUE) +
    log(truncnorm::dtruncnorm(a.new, a = lb, b = ub, mean = (lb+ub)/2, sd = 2)) -
    log(truncnorm::dtruncnorm(a.old, a = lb, b = ub, mean = (lb+ub)/2, sd = 2)) +
    log(d.old) - log(d.new);
  }

  if (runif(1) < exp(log.accept)) alpha = a.new

  return(alpha)
}

## this is for alpha = (0,1) if alpha>1, dbeta(a.new, pr.a, pr.b, log=TRUE) = -Inf and log.accept = NaN
draw.alpha.randomwalk <- function(alpha, beta, tau, window, pr.a=1.0, pr.b=1.0, ep=0.1)
{
                                        # if (sum(abs(beta) <= 1e-5) == 0 & abs(tau) > 1e-5) {
    s = log(abs(beta / tau));
                                        # } else {
                                        #   beta[abs(beta) <= 1e-5] <- 1e-4
                                        #   if (abs(tau) <= 1e-5) tau <- 1e-4
                                        #   s = log(abs(beta / tau));
                                        # }

    a.old = alpha;
    lb = 0.0
    ub = 2.0

    a.new = a.old + rnorm(1, 0, sd=sqrt(window));
    if(a.new > 2){
        a.new = min(c(ub, a.new+ep))
    }else if(a.new < 0){
        a.new = max(c(lb, a.new-ep))
    }else{
        a.new = a.new
    }
    ## d.old = min(c(ub, a.new+ep)) -  max(c(lb, a.new-ep))

    log.accept = llh.alpha(a.new, s) - llh.alpha(a.old, s)
    ## +  dbeta(a.new, pr.a, pr.b, log=TRUE) - dbeta(a.old, pr.a, pr.b, log=TRUE);

    if (runif(1) < exp(log.accept)) {
        alpha = a.new
        cat("\n accepted alpha is ", alpha, "\n")
    }
    return(alpha)
}




## this is for alpha = (0,1) if alpha>1, dbeta(a.new, pr.a, pr.b, log=TRUE) = -Inf and log.accept = NaN
draw.alpha2 <- function(alpha, beta, tau, pr.a=1.0, pr.b=1.0, ep=0.1)
{
  s = log(abs(beta / tau));
  a.old = alpha;
  lb = 0.0
  ub = 1.0

  ## We might want to change these bounds if we have some belief about the value
  ## of alpha.
  l.new = max(c(lb, a.old-ep))
  r.new = min(c(ub, a.old+ep))
  d.new = r.new - l.new
  a.new = runif(1, l.new, r.new);

  d.old = min(c(ub, a.new+ep)) -  max(c(lb, a.new-ep))

  ## s.new = tau * a.new^(-1/a.new);
  ## s.old = tau * a.old^(-1/a.old);
  ##log.accept = sum(mydpgnorm(beta, 0, tau, a.new, log=TRUE)) - sum(mydpgnorm(beta, 0, tau, a.old, log=TRUE)) +
    log.accept = llh.alpha(a.new, s) - llh.alpha(a.old, s) +
    dbeta(a.new, pr.a, pr.b, log=TRUE) - dbeta(a.old, pr.a, pr.b, log=TRUE) +
    log(d.old) - log(d.new);

  if (runif(1) < exp(log.accept)) alpha = a.new

  alpha
}



dtruncnorm <- function(x, lb, ub, log = TRUE) {
    mean <- (ub + lb) / 2
    if (log) {
        density <- dnorm(x, mean = mean, sd = sqrt(3), log = TRUE) -
                    log(pnorm(ub, mean = mean, sd = sqrt(3)) - pnorm(lb, mean = mean, sd = sqrt(3)))
    } else {
        density <- exp(dnorm(x - mean, log = TRUE) - log(pnorm(ub - mean) - pnorm(lb - mean)))
    }
    return(density)
}

draw.alpha.griddy <- function(beta, tau, lb = 0.0, ub = 2.0) {

    s = log(abs(beta / tau))
    # lb = 0.0; ub = 2.0 ## upper bound must be at least 2 (ridge estimator).

    grid <- seq(lb+0.1, ub, by = 0.01)
    density <- rep(NA, length(grid))
    for (j in 1:length(grid)) {
      # density[j] <- llh.alpha(grid[j], s) + dbeta(grid[j], pr.a, pr.b, log = TRUE)
      # density[j] <- llh.alpha(grid[j], s) + dunif(grid[j], lb, ub, log = TRUE)
        density[j] <- llh.alpha(grid[j], s) +
            log(truncnorm::dtruncnorm(grid[j], lb, ub, mean = (lb+ub)/2, sd = 1.5))
    }

    ## if (sum(exp(density) == 0) == length(density)){
                                        # require(Rmpfr)
    ## density <- as(density, "mpfr")
    ## detach("package:Rmpfr", unload=TRUE)
    ##}

    ## slice
    tryCatch(alpha <- sample(grid, size = 1, prob = exp(density)),
             warning = function(w) {cat("Warning: alpha griddy sampler density is ", exp(density))},
             error = function(e) {cat("Error: alpha griddy sampler density is ", exp(density))}
              )
    ## Error in sample.int(length(x), size, replace, prob) :
    return(alpha)
}



## from BayesLogit
draw.df <- function(d.prev, psi, G, ymax)
{
  d.new = d.prev

  ## Kernel 1: RW MH.
  ## And a random walk.
  rw.lower = max(d.prev - 1, 1);
  rw.upper = d.prev + 1;
  rw.grid  = rw.lower:rw.upper;
  rw.n     = length(rw.grid)
  rw.p     = rep(1/rw.n, rw.n);

  d.prop = sample(rw.grid, 1, prob=rw.p);

  ltarget = df.llh(d.prop, psi, G, ymax) - df.llh(d.prev, psi, G, ymax)
  lpps = log(ifelse(d.prop==1, 1/2, 1/3)) - log(ifelse(d.prev==1, 1/2, 1/3));
  lratio = ltarget + lpps

  if (runif(1) < exp(lratio)) d.new = d.prop

  d.new
}

df.llh <- function(d, psi, G, ymax)
{
  p =  1 / (1 + exp(-psi))
  sum( log(d+0:(ymax-1)) * G[1:ymax] ) + d * sum(log(1-p));
}


draw.tau <- function(beta, alpha, c, d)
{
  p = length(beta)
  nu = rgamma(1, c + p/alpha, rate=d + sum(abs(beta)^alpha))
  tau = nu^(-1/alpha)
  return(tau);
}

draw.sig2 <- function(beta, x, y, sig2.shape, sig2.scale)
{
  n = length(y)
  rss = sum( (as.matrix(y)-x%*%as.matrix(beta))^2 )
  prec = rgamma(1, shape = sig2.shape+n/2, rate = sig2.scale+rss/2)
  return(1/prec)
}
## betaj = beta[j,]; bij= bi[j,]; ## , Xj, yj, Wj,  sig2.shape=c0, sig2.scale=d0
draw.panel.sig2 <- function(betaj, bij, Xj, yj, Wj, sig2.shape=0.0, sig2.scale=0.0)
{
    nj = length(yj)
    rss = t(as.matrix(yj)-Xj%*%as.matrix(betaj)-Wj%*%as.matrix(bij))%*%
        (as.matrix(yj)-Xj%*%as.matrix(betaj)-Wj%*%as.matrix(bij))
    ## sum( (as.matrix(y)-x%*%as.matrix(beta))^2 )
    prec = rgamma(1, sig2.shape+nj/2, rate=sig2.scale+rss/2)
    return(1/prec)
}

sparse.state.sampler <- function(m=m, y=y, X=X, beta=beta, sig2=sig2, P = P){
        T   <-  length(y)
        F   <-  matrix(NA, T, m+1)     # storage for the Filtered probabilities
        pr1 <-  c(1,rep(0, m))         # initial probability Pr(s=k|Y0, lambda)
        for (t in 1:T){
            mu <- X[t,] %*% t(beta)
            py  <-  sapply(c(1:(m+1)), function(i){
                dnorm(y[t], mean= mu[i], sd = sqrt(sig2[i]))})
            if(t==1) {
                pstyt1 = pr1
            } else {
                pstyt1     <- F[t-1,]%*%P
            }
            unnorm.pstyt    <- pstyt1*py
            pstyt   <-  unnorm.pstyt/sum(unnorm.pstyt) # Pr(st|Yt)
            F[t,]   <-  pstyt
        }

        s      <-  matrix(NA, T, 1)   ## holder for state variables
        ps     <-  matrix(NA, T, m+1) ## holder for state probabilities
        ps[T,] <-  F[T,]              ## we know last elements of ps and s
        s[T,1] <-  m+1
        s[1] <- 1
        ## t      <-  T-1
        for(t in (T-1):2){
            ## while (t>=1){
            st     <-  s[t+1]
            unnorm.pstyn   <-  F[t,]*P[,st]
            if(sum(unnorm.pstyn) == 0){
                ## if unnorm.pstyn is all zero, what to do?
                cat("state sampler stops at t = ", t, " because F", F[t,]," and P", P[,st]," do not match! \n")
                s[t]   <- s[t+1]
            }else{
                ## normalize into a prob. density
                pstyn   <-  unnorm.pstyn/sum(unnorm.pstyn)
                if (st==1) {
                    s[t]<-1
                } else {
                    pone    <-  pstyn[st-1]
                    s[t]   <-  ifelse (runif(1) < pone, st-1, st)
                }
                ps[t,] <-  pstyn
                ## probabilities pertaining to a certain state
            }
        }

        ## name and report outputs
        out <-  list(s, ps)
        names(out) <- c("state","ps")
        return(out)
    }

sparse.state.sampler.NB <- function(m=m, y=y, d=d, X=X, beta=beta, P = P){
        T   <-  length(y)
        F   <-  matrix(NA, T, m+1)     # storage for the Filtered probabilities
        pr1 <-  c(1,rep(0, m))         # initial probability Pr(s=k|Y0, lambda)

        for (t in 1:T){
            ## for NB
            exb <- exp(X[t,] %*% t(beta))
            prob1 <- exb/(1+exb)
            py  <-  sapply(c(1:(m+1)), function(i){
                dbinom(y[t], size = y[t] + d, prob = prob1[i])})

            if(t==1) {
                pstyt1 = pr1
            } else {
                pstyt1     <- F[t-1,]%*%P
            }
            unnorm.pstyt    <- pstyt1*py
            pstyt   <-  unnorm.pstyt/sum(unnorm.pstyt) # Pr(st|Yt)
            F[t,]   <-  pstyt
        }

        s      <-  matrix(NA, T, 1)   ## holder for state variables
        ps     <-  matrix(NA, T, m+1) ## holder for state probabilities
        ps[T,] <-  F[T,]              ## we know last elements of ps and s
        s[T,1] <-  m+1
        s[1] <- 1
        ## t      <-  T-1
        for(t in (T-1):2){
            ## while (t>=1){
            st     <-  s[t+1]
            unnorm.pstyn   <-  F[t,]*P[,st]
            if(sum(unnorm.pstyn) == 0){
                ## if unnorm.pstyn is all zero, what to do?
                cat("state sampler stops at t = ", t, " because F", F[t,]," and P", P[,st]," do not match! \n")
                s[t]   <- s[t+1]
            }else{
                ## normalize into a prob. density
                pstyn   <-  unnorm.pstyn/sum(unnorm.pstyn)
                if (st==1) {
                    s[t]<-1
                } else {
                    pone    <-  pstyn[st-1]
                    s[t]   <-  ifelse (runif(1) < pone, st-1, st)
                }
                ps[t,] <-  pstyn
                ## probabilities pertaining to a certain state
            }
        }

        ## name and report outputs
        out <-  list(s, ps)
        names(out) <- c("state","ps")
        return(out)
    }
## centering data
## m=m, y=y, X=X, W=W, beta=beta, bi=bi, sig2=sig2, D=D, P=P)

sparse.panel.state.sampler.NB <- function(m=m, d=d, Yt_arr=Yt_arr, Xt_arr=Xt_arr, Wt_arr=Wt_arr,
                                          beta=beta, bi=bi, P = P){
    ## T   <-  length(y)
    F   <-  matrix(NA, T, m+1)     # storage for the Filtered probabilities
    pr1 <-  c(1,rep(0, m))         # initial probability Pr(s=k|Y0, lambda)
    py <- rep(NA, m+1)
    for (tt in 1:T){
        yj <- Yt_arr[[tt]]
        for(j in 1:(m+1)){
            ## for panel NB
            mu.state  <-  Xt_arr[[tt]]%*%beta[j,] + Wt_arr[[tt]]%*%bi[[j]]
            exb  <-  exp(mu.state)
            prob1 <- exb/(1+exb)
            py[j]  <-  prod(dbinom(yj, size = yj + d, prob = prob1))
        }
        if(tt==1) {
            pstyt1 = pr1
        }else {
            pstyt1     <- F[tt-1,]%*%P
        }
        unnorm.pstyt    <- pstyt1*py
        pstyt   <-  unnorm.pstyt/sum(unnorm.pstyt) # Pr(st|Yt)
        F[tt,]   <-  pstyt
    }

    s      <-  matrix(NA, T, 1)   ## holder for state variables
    ps     <-  matrix(NA, T, m+1) ## holder for state probabilities
    ps[T,] <-  F[T,]              ## we know last elements of ps and s
    s[T,1] <-  m+1
    s[1] <- 1
    ## t      <-  T-1
    for(t in (T-1):2){
        ## while (t>=1){
        st     <-  s[t+1]
        unnorm.pstyn   <-  F[t,]*P[,st]
        if(sum(unnorm.pstyn) == 0){
            ## if unnorm.pstyn is all zero, what to do?
            cat("state sampler stops at t = ", t, " because F", F[t,]," and P", P[,st]," do not match! \n")
            s[t]   <- s[t+1]
        } else{
            ## normalize into a prob. density
            pstyn   <-  unnorm.pstyn/sum(unnorm.pstyn)
            if (st==1) {
                s[t]<-1
            }else {
                pone    <-  pstyn[st-1]
                s[t]   <-  ifelse (runif(1) < pone, st-1, st)
            }
            ps[t,] <-  pstyn
            ## probabilities pertaining to a certain state
        }
    }

    ## name and report outputs
    out <-  list(s, ps)
    names(out) <- c("state","ps")
    return(out)
}


## --> replaced with sparse.panel.state.sampler.cpp
sparse.panel.state.sampler <- function(m=m, T=T, N=N, Yt_arr=Yt_arr, Xt_arr=Xt_arr,
                                       Wt_arr=Wt_arr, beta=beta, bi=bi, sig2=sig2, D=D, P=P){
    ## T   <-  length(y)
    F   <-  matrix(NA, T, m+1)     # storage for the Filtered probabilities
    pr1 <-  c(1,rep(0, m))         # initial probability Pr(s=k|Y0, lambda)
    py <- rep(NA, m+1)
    for (tt in 1:T){
        for(j in 1:(m+1)){
            Mu <- matrix(Xt_arr[[tt]] %*% beta[j, ], N, 1)
            Sigma <- sig2[j]*diag(N) + Wt_arr[[tt]]%*%D[[j]]%*%t(Wt_arr[[tt]])
            py[j] <- dmvnorm(t(Yt_arr[[tt]]), mean= Mu, sigma = Sigma)
        }
        if(tt==1) {
            pstyt1 = pr1
        }else {
            pstyt1     <- F[tt-1,]%*%P
        }
        unnorm.pstyt    <- pstyt1*py
        pstyt   <-  unnorm.pstyt/sum(unnorm.pstyt) # Pr(st|Yt)
        F[tt,]   <-  pstyt
    }

    s      <-  matrix(NA, T, 1)   ## holder for state variables
    ps     <-  matrix(NA, T, m+1) ## holder for state probabilities
    ps[T,] <-  F[T,]              ## we know last elements of ps and s
    s[T,1] <-  m+1
    s[1] <- 1
    ## t      <-  T-1
    for(t in (T-1):2){
        ## while (t>=1){
        st     <-  s[t+1]
        unnorm.pstyn   <-  F[t,]*P[,st]
        if(sum(unnorm.pstyn) == 0){
            ## if unnorm.pstyn is all zero, what to do?
            cat("state sampler stops at t = ", t, " because F", F[t,]," and P", P[,st]," do not match! \n")
            s[t]   <- s[t+1]
        } else{
            ## normalize into a prob. density
            pstyn   <-  unnorm.pstyn/sum(unnorm.pstyn)
            if (st==1) {
                s[t]<-1
            }else {
                pone    <-  pstyn[st-1]
                s[t]   <-  ifelse (runif(1) < pone, st-1, st)
            }
            ps[t,] <-  pstyn
            ## probabilities pertaining to a certain state
        }
    }

    ## name and report outputs
    out <-  list(s, ps)
    names(out) <- c("state","ps")
    return(out)
}
sparse.panel.state.sampler2 <- function(m=m, T=T, N=N, Yt_arr=Yt_arr, Xt_arr=Xt_arr,
                                       beta=beta, br=br, sig2=sig2, D=D, P=P){
    ## T   <-  length(y)
    F   <-  matrix(NA, T, m+1)     # storage for the Filtered probabilities
    pr1 <-  c(1,rep(0, m))         # initial probability Pr(s=k|Y0, lambda)
    py <- rep(NA, m+1)
    for (tt in 1:T){
        for(j in 1:(m+1)){
            Mu <- matrix(Xt_arr[[tt]] %*% beta[j, ] + t(br[[j]]), N, 1)
            Sigma <- sig2[j]*diag(N) ## + Wt_arr[[tt]]%*%D[[j]]%*%t(Wt_arr[[tt]])
            py[j] <- dmvnorm(t(Yt_arr[[tt]]), mean= Mu, sigma = Sigma)
        }
        if(tt==1) {
            pstyt1 = pr1
        }else {
            pstyt1     <- F[tt-1,]%*%P
        }
        unnorm.pstyt    <- pstyt1*py
        pstyt   <-  unnorm.pstyt/sum(unnorm.pstyt) # Pr(st|Yt)
        F[tt,]   <-  pstyt
    }

    s      <-  matrix(NA, T, 1)   ## holder for state variables
    ps     <-  matrix(NA, T, m+1) ## holder for state probabilities
    ps[T,] <-  F[T,]              ## we know last elements of ps and s
    s[T,1] <-  m+1
    s[1] <- 1
    ## t      <-  T-1
    for(t in (T-1):2){
        ## while (t>=1){
        st     <-  s[t+1]
        unnorm.pstyn   <-  F[t,]*P[,st]
        if(sum(unnorm.pstyn) == 0){
            ## if unnorm.pstyn is all zero, what to do?
            cat("state sampler stops at t = ", t, " because F", F[t,]," and P", P[,st]," do not match! \n")
            s[t]   <- s[t+1]
        } else{
            ## normalize into a prob. density
            pstyn   <-  unnorm.pstyn/sum(unnorm.pstyn)
            if (st==1) {
                s[t]<-1
            }else {
                pone    <-  pstyn[st-1]
                s[t]   <-  ifelse (runif(1) < pone, st-1, st)
            }
            ps[t,] <-  pstyn
            ## probabilities pertaining to a certain state
        }
    }

    ## name and report outputs
    out <-  list(s, ps)
    names(out) <- c("state","ps")
    return(out)
}

draw.u <- function(tau, beta, w, alpha)
{
  m = 1-{abs(beta/tau)*w^(-1/alpha)}
  runif(length(beta),max=m)
}

draw.w <- function(tau, beta, u, alpha)
{
  p = length(beta)
  a = (abs(beta/tau)/(1-u))^alpha
  pr = alpha/(1+alpha*a)
  ## p1 = 0.5*(1+alpha)
  ## pr = p1/(1-p1*a)
  shape = (runif(p) < pr) + 1;
  w = rgamma(p, shape, 1);
  w = w+a
  list("w"=w, "shape"=shape);
}

draw.w.marg <- function(alpha, p)
{
  pr = 1-alpha
  shape = (runif(p) < pr) + 1;
  w = rgamma(p, shape, 1)
  list("w"=w, "shape"=shape);
}


##------------------------------------------------------------------------------

sig.for.pg <- function(tau, alpha)
{
  sig2 = tau^2 * gamma(3/alpha) / gamma(1/alpha);
  sqrt(sig2)
}

## mydpgnorm <- function(y, m, tau, alpha, log=FALSE)
## {
##   sig = sig.for.pg(tau, alpha);
##   out = dpgnorm(y, alpha, m, sig);
##   if (log) out = log(out)
##   out
## }


draw.alpha.2 <- function(alpha, beta, tau, ep=0.1)
{
  s = log(abs(beta / tau));
  a.old = alpha;
  lb = 0.01
  ub = 0.99

  ## We might want to change these bounds if we have some belief about the value
  ## of alpha.
  a.new = runif(1, lb, ub);

  d.old = min(c(ub, a.new+ep)) -  max(c(lb, a.new-ep))

  log.accept = - sum(exp(s*a.new)) + sum(exp(s*a.old))

  if (runif(1) < exp(log.accept)) alpha = a.new

  alpha
}

##------------------------------------------------------------------------------
## Sampling beta

## Geweke style Gibbs sampling
draw.beta.1 <- function(beta, bhat, xx, sig2, tau, u, w, alpha)
{
  p = length(bhat)
  b = (1-u)*{w^(1/alpha)}*tau
  for ( i in 1:p )
    {
      m = bhat[i] - crossprod(xx[i,-i],beta[-i]-bhat[-i])/xx[i,i]
      v = sig2/xx[i,i]
      beta[i] = rtnorm(1,m,sqrt(v),-b[i],b[i])
    }
  beta
}

## Rodriguez-Yam style Gibbs sampling - using SVD
draw.beta.2 <- function(beta, a, tV, d, sig2, tau, u, w, alpha)
{
  P = length(beta);

  b = (1-u)*{w^(1/alpha)}*tau
  ## b = (1-u) * tau * exp(log(w) / alpha)
  z = tV %*% beta;

  for (i in 1:P) {
    lmax = -Inf
    rmin =  Inf

    for (j in 1:P) {
      vji = tV[i,j];
      vj  = tV[ ,j];
      ## rji = vj %*% z - vji * z[i];
      rji = vj[-i] %*% z[-i];
      Dif = b[j] - rji
      Sum = b[j] + rji
      left  = ifelse(vji > 0, -Sum, -Dif) / abs(vji);
      right = ifelse(vji > 0,  Dif,  Sum) / abs(vji);
      lmax = max(c(lmax,left))
      rmin = min(c(rmin,right))
    }

    if (d[i]!=0) {
      m = a[i] / (d[i]^2);
      s = sqrt(sig2) / d[i];
      z[i] = rtnorm(1, m, s, lmax, rmin);
    } else {
      cat("d =", d[i], "\n");
      z[i] = runif(1, lmax, rmin);
    }

  }

  beta = t(tV) %*% z;
  beta
}

## Alternate Rodriguez-Yam style Gibbs sampling.  Not as good.
draw.beta.3 <- function(beta, bhat, xx, sig2, tau, u, w, alpha)
{
  P = length(beta)
  evd = eigen(xx)
  rt  = evd$vectors %*% diag(sqrt(evd$values), P) %*% t(evd$vectors);  ## sqrt XX
  irt = evd$vectors %*% diag(sqrt(1/evd$values), P) %*% t(evd$vectors);
  b = (1-u)*{w^(1/alpha)}*tau

  z = rt %*% beta
  m = rt %*% bhat

  for (i in 1:P) {

    left  = rep(0, P)
    right = rep(0, P)

    for (j in 1:P) {
      rji = irt[j,-i] %*% z[-i]
      Dif = b[j] - rji;
      Sum = b[j] + rji;
      left[j]  = ifelse(irt[j,i] > 0, -Sum, -Dif) / abs(irt[j,i]);
      right[j] = ifelse(irt[j,i] > 0,  Dif,  Sum) / abs(irt[j,i]);
    }

    lmax = max(left)
    rmin = min(right)

    z[i] = rtnorm(1, m[i], sqrt(sig2), lmax, rmin);

  }

  beta = irt %*% z;
}



##
## SLOG function from Rajaratnam et al (2016; JRSSB)
## this function is used to initialize Beta's
##
SLOG <- function(x, y, l, times = 1e-6, thresh = 1e-10, start=NULL){
    #function for implementing SLOG/rSLOG
    #x: covariate data
    #y: response data
    #l: value of lasso regularizaton parameter
    #times: convergence criteria - difference between successive coefficient vectors
    #thresh: below this value estimates are set to 0 (runs rSLOG)
    #start: allows starting values other than sign(xty)*l/p to be specified.
    xtx <- crossprod(x)
    xty <- crossprod(x,y)
    p <- length(xty)
    n <- length(y)
    b.cur <- sign(xty) * l / p
    if(!is.null(start)) b.cur <- start
        b.old <- b.cur
        vin  <- 1:p
        temp <-rep(0,p)
        conv = FALSE
        k <- 1
    while(conv==FALSE){
        b <- b.cur[vin]
        p <- length(b)
        B.inv <- diag(l/abs(b),nrow=p)
        b.new <- tcrossprod(chol2inv(chol(B.inv+xtx[vin,vin])),t(xty[vin]))
        temp[vin] <- as.vector(b.new)
        # temp[abs(temp) <= thresh] <- 0
        b.new <- temp
        vin   <- which(b.new!=0)
        b.cur <- as.vector(b.new)
        conv  <- (sqrt(sum((b.cur-b.old)^2))/sqrt(sum(b.old^2)))<times
        b.old <- b.cur
        k <- k+1
    }

    return(b.cur)
}
## ## Hamiltonian MCMC
## draw.beta.4 <- function(beta, bhat, xx, sig2, tau, u, w, alpha)
## {
##   require("tmg")
##   p = length(beta);

##   b = (1-u)*{w^(1/alpha)}*tau

##   ## Constraints
##   F = matrix(nrow=2*p, ncol=p)
##   F[1:p,1:p]   = diag(1,p)
##   F[1:p+p,1:p] = -1 * diag(1,p)
##   g = as.vector(c(b,b))

##   ## print(b)
##   ## print(beta)

##   prec = xx / sig2;
##   ell  = prec %*% bhat;

##   out = rtmg(1, prec, ell, beta, F, g, burn.in = 30);

##   drop(out)
## }


#' Function to update intercept 
#' @param y outcome vector 
#' @param Xorig design matrix in original scale 
#' @param beta regression coefficient
#' @param n.break number of breaks 
#' @param intercept boolean; if \code{FALSE}, the intercept is set to zero.
#' @param state state vector (length T)
#' @keywords internal
estimate_intercept_reg <- function(y, Xorig, beta, n.break, intercept, state) {
  ns <- n.break + 1

  if (n.break == 0 & intercept == TRUE) {
    ## fit intercept
    beta0 <- mean(y) - as.vector(colMeans(Xorig) %*% t(beta))
    # cat("beta0 = ", beta0, "\n")
  } else if (n.break > 0 & intercept == TRUE) {
    beta0 <- matrix(NA, nrow = ns, ncol = 1)
    # ydm   <- as.vector(as.vector(y) - tapply(y, state, mean)[state])
    for (j in 1:ns) {
      beta0[j,] <- mean(y[state == j]) - colMeans(Xorig[state == j, , drop = FALSE]) %*% beta[j,]
    }
  } else if (intercept == FALSE) {
    beta0 <- matrix(0, nrow = ns, ncol = 1)
    # ydm   <- as.vector(as.vector(y) - tapply(y, state, mean)[state])
  }

  return(beta0)
}
