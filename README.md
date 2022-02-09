
<!-- README.md is generated from README.Rmd. Please edit that file -->
BridgeChange
============

`R` package `BridgeChange` fits hidden Markove Bayesian bridge models for time-series data and panel data with possibly a large number of covariates and change-points. This package can be used to detect, estimate, and summarize time-varying parameters of regression models where high-dimensinal covariates have heterogenous effects on the response data across time. The fixed-effects model, the random-effects model (random coefficient model or multilevel model), and the univariate regression model are available. Model diagnostics can be done by Watanabe-Akaike Information Criterion (WAIC). Posterior estimates can be summarized by the decoupled shrinkage and selection (DSS) method. 

Table of Contents
-----------------

1.  [Overview](#overview)
2.  [Installation](#installation)
3.  [Time-series Data](#time-series-data)

Installation
------------

You can install the most recent version of `BridgeChange` from Gitub using the [`devtools`](https://github.com/r-lib/devtools) package.

``` r
# install BridgeChange from Github
# you might need to instal "devtools"
devtools::install_github("jongheepark/BridgeChange")
```

Time-series Data
----------------

### Model

`BridgeChangeReg()` can be used to analyse time-series data with possible change-points.

The functions fits the linear model:
*y*<sub>*t*</sub> = **X**<sub>*t*</sub><sup>⊤</sup>*β*<sub>*s*<sub>*t*</sub></sub> + *ϵ*<sub>*s*<sub>*t*</sub></sub>
 where *s*<sub>*t*</sub> ∈ {1, …, *M*} is an indicator of states.

### Example

``` r
set.seed(1973);

library("pcse")
set.seed(1999)
data("agl")
data = agl
model = "within"
index = c('country', 'year')
effect = 'time'
formula <- growth ~ lagg1 + opengdp + openex + openimp + leftc + central + inter
```

``` r
pdata   <- pdata.frame(data, index)
pm <- plm(formula, data = pdata, model = model, effect = effect)

## plot of panel residuals 1
coplot(pm$residuals ~ pdata[,index[2]]|pdata[,index[1]], data=pdata, 
       overlap=.1, col="brown", type="l", panel = panel.smooth, xlab="panel residuals by group and time")
```


``` r
##
## fit models 
##

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

```

``` r
##
## Post-estimation 
##

## model selection by WAIC
waic <- WaicCompare(list(agl.cp0, agl.cp1, agl.cp2), print = TRUE)
plotWaic(waic)

par(mfrow=c(1, 2))
MCMCpack::plotState(agl.cp1, start=1970, legend.control =c(1970, 0.85), main="One break")
MCMCpack::plotState(agl.cp2, start=1970, legend.control =c(1970, 0.85), main="Two breaks")
```
