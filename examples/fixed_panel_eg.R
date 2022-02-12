####################################
## Alvarez et al. data
####################################
library("pcse")
set.seed(1999)
data("agl")
data = agl
model = "within"
index = c('country', 'year')
effect = 'time'
formula <- growth ~ lagg1 + opengdp + openex + openimp + leftc + central + inter



#########################
## inspection of residuals
#########################
pdata   <- pdata.frame(data, index)
pm <- plm(formula, data = pdata, model = model, effect = effect)
                                        # 
## plot of panel residuals 1
coplot(pm$residuals ~ pdata[,index[2]]|pdata[,index[1]], data=pdata, ## number=length(unique(pdata[,index[1]])),
       overlap=.1, col="brown", type="l", panel = panel.smooth, xlab="panel residuals by group and time")

## plot of panel residuals 2
coplot(pm$residuals ~ pdata[,index[2]]|pdata[,index[1]], type="b", data=pdata)


#########################
## HMBB estimation
#########################
## Be aware that 100 runs are too short for analysis. 
mcmc = 100; burn = 100; verbose = 100; thin = 1;
formula <- growth ~ lagg1 + opengdp + openex + openimp + leftc + central + inter
agl.cp0 <- BridgeFixedPanel(formula=formula, data = data, model = model, index = index, effect = effect,
                            mcmc=mcmc, verbose=verbose, Waic = TRUE, 
                            n.break = 0)
agl.cp1 <- BridgeFixedPanel(formula=formula, data = data, model = model, index = index, effect = effect,
                            mcmc=mcmc, verbose=verbose, Waic = TRUE, 
                            n.break = 1)
agl.cp2 <- BridgeFixedPanel(formula=formula, data = data, model = model, index = index, effect = effect,
                            mcmc=mcmc, verbose=verbose, Waic = TRUE, 
                            n.break = 2)

## model selection by WAIC
waic <- WaicCompare(list(agl.cp0, agl.cp1, agl.cp2), print = TRUE)
plotWaic(waic)

par(mfrow=c(1, 2))
MCMCpack::plotState(agl.cp1, start=1970, legend.control =c(1970, 0.85), main="One break")
MCMCpack::plotState(agl.cp2, start=1970, legend.control =c(1970, 0.85), main="Two breaks")

####################################
## dotplot over time
## time-varying movements of selected covariates
####################################
## all covariates
dotplotRegime(agl.cp1, hybrid=FALSE, start = 1970, location.bar=12, x.location="default",
              text.cex=0.8, main="Time-varying Movements of All Covariates")
## label as a legend
dotplotRegime(agl.cp1, hybrid=FALSE, start = 1970, location.bar=12, x.location="legend",
              text.cex=0.8, main="Time-varying Movements of All Covariates")

## leftc only
## select works like grep()
dotplotRegime(agl.cp1, hybrid=FALSE, start = 1970, location.bar=12, x.location="static",
              text.cex=0.8, select="left", main=("Time-varying Movements of Left Party"))
