
## hybrid
## location.bar is the x location of the legend starting
## select = grep coefficients containing a string
## hybrid=TRUE and raw=TRUE, use original y for variable selection
## hybrid=TRUE and raw=FALSE, use yhat for variable selection
## order = k : Order the transparency of coefficients based on regime k. The default is k=1.
## smooth.beta If TRUE, beta*pr(state==j) is used to compute time beta for state j.
##              Otherwise, beta*I(state==j) is used where I() is an indicator function. 
## dotplotRegime <- function(out, hybrid=TRUE, start=1, cex=1, smooth.beta = TRUE, 
##                           x.location=c("random", "legend", "default"),
##                           order.state = 1, location.bar=9,
##                           text.cex=1, legend.position = "topright", pos=1, ## below default, 3 above
##                          select=NULL, main="", raw=FALSE){

#' Draw a plot of covariate-specific time-varying movements 
#'
#'
#' @param out output of HMBB
#' @param hybrid If TRUE, DSS results are used for graph. Otherwise, HMBB results are used. 
#' @param start Starting value of the time sequence
#'
#' @param smooth.beta If TRUE, probabilistic estimates of hidden states (smooth) are used. Otherwise, discrete state values are used. 
#'
#' @param x.location x location of covariate name labels. If random, a random value (between min and max of x) is selected.
#' If legend, it will be displayed in the top-right as a legend. If default, the labels will be shown in the above of the last regime estimates. 
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
#' @param h hue value in the HCL color description, has to be in [0, 360]. Default is c(255, 330).
#'
#' @param l luminance value in the HCL color description. Default is c(40, 90).
#'
#' @param distance Distance of text labels from the maximum value of x. Default is 4. 
#' 
#' @export
#'

dotplotRegime <- function(out, hybrid=TRUE, start=1, cex=1, smooth.beta = TRUE, 
                          x.location=c("none", "random", "legend", "default"),
                          order.state = 1, location.bar=9, distance = 4, 
                          h=c(255, 330), l = c(40, 90),
                          text.cex=1, legend.position = "topright",
                          pos=1, ## below default, 3 above
                          select=NULL, main="", raw=FALSE){
    if(attr(out, "m") == 0){
        state <- rep(1, length(attr(out, "y")))
    }else{
        state <- round(apply(attr(out, "s.store"), 2, mean))
    }
    unique.time.index <- start : (start + length(state) - 1)
    m <- attr(out, "m")
    ns <- m + 1
    X <- attr(out, "X")
    k <- dim(X)[2]
    if(hybrid){
        if(raw){
            coef <- attr(out, "hybrid.raw")
        }else{
            coef <- attr(out, "hybrid")
        }
    }else{
        ### coefs <- summary(out)[[1]][,1]
        beta.target <- out[, grepl("beta", colnames(out))]
        coef <- matrix(apply(beta.target, 2, mean), k, ns)
        rownames(coef) <- colnames(X)     
    }
    ## select
    if(!is.null(select)){
        if(length(grep(select, rownames(coef))) == 1){
            coefs <-  matrix(coef[ grep(select, rownames(coef)), ], 1, m+1)
            rownames(coefs) <- rownames(coef)[grep(select, rownames(coef))]
        }else{
            coefs <- coef[grep(select, rownames(coef)),]
        }
    }else{
        coefs <- coef
    }


    ## if both regime estimates are zero, drop them.
    if(sum(apply(coef, 1, prod) + apply(coef, 1, sum) == 0)>0){
        coef <- coef[-which(apply(coef, 1, prod) + apply(coef, 1, sum) == 0),]
    }
    coef.mat <- matrix(NA, nrow=nrow(coefs), ncol=length(unique.time.index))
    
    if(smooth.beta){
        prob.state <- cbind(sapply(1:ns, function(k){apply(attr(out, "s.store") == k, 2, mean)}))
        coef.mat <- coefs%*%t(prob.state)
    }else{
        for(i in 1:nrow(coefs)){
            coef.mat[i,] <- coefs[i, state]
        }
    }

    
    ## require(colorspace)
    if(nrow(coef.mat) == 1){
        col.scheme <- "brown"
    }else{
        col.scheme <- diverge_hcl(nrow(coefs), h=h, l =l)
        col.scheme <- col.scheme[rank(coefs[, order.state])]
    }
    ## par (mar=c(3,3,2,4), mgp=c(2,.7,0), tck=-.01)
    plot(unique.time.index, coef.mat[1,], xlab="time", ylab="coefficients", bty='n', main=main,
         ylim=range(coef.mat), type="o", pch=19, cex=cex, col=addTrans(col.scheme[1], 100))
    if(nrow(coef.mat)>1){
        for(i in 2:nrow(coefs)){
            points(unique.time.index, coef.mat[i,], pch=19, cex=cex, col=addTrans(col.scheme[i], 100))
            lines(unique.time.index, coef.mat[i,], lwd=0.8, col=col.scheme[i])
        }
    }
    ## grid(col="grey80")
    if(x.location=="random"){
        n.coef <- length(rownames(coefs))
        if(min(unique.time.index)+location.bar != (max(unique.time.index)-2)){
            x.location.pos <- sample((min(unique.time.index)+location.bar):(max(unique.time.index)-2),
                                     n.coef, replace=TRUE)
        }else{
            x.location.pos <- min(unique.time.index)+location.bar
        }
        text(x = x.location.pos, coef.mat[, length(unique.time.index)], rownames(coefs),
             cex=text.cex, col=col.scheme,
             pos=pos)
    }else if(x.location=="default"){
        text(x = max(unique.time.index)- distance, coef.mat[, length(unique.time.index)], pos=pos, 
             rownames(coefs), col=col.scheme, cex=text.cex)

    }else if(x.location=="legend"){
        legend(legend.position, legend=rownames(coefs), col=addTrans(col.scheme, 100), pch=19, bty="n",
               cex=0.8,lty=1:1, lwd=1, y.intersp = 0.8)
    }else{
    }
    ## title("Regime-changing coefficients", adj = 0, line = 0)
    box(); 
}
