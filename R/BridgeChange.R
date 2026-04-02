#' @useDynLib BridgeChange, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom MCMCpack rinvgamma rwish dinvgamma dwish
#' @importFrom mvtnorm rmvnorm dmvnorm
#' @importFrom copula retstable
#' @importFrom coda as.mcmc mcmc
#' @importFrom plm plm pmodel.response index
#' @importFrom colorspace diverge_hcl
#' @importFrom truncnorm dtruncnorm
#' @importFrom glmnet cv.glmnet
#' @importFrom ggplot2 ggplot aes geom_line geom_point facet_wrap labs theme_bw theme element_blank element_text
#' @importFrom grDevices col2rgb
#' @importFrom graphics abline axis box legend lines par points polygon text
#' @importFrom stats as.formula ave cor dbeta dbinom dgamma dnorm fitted median pnorm qnorm quantile rexp rgamma rnorm runif ts vcov coef complete.cases lm model.matrix sd var
#' @importFrom utils combn flush.console getFromNamespace
NULL
#> NULL

## Suppress R CMD check NOTEs for non-standard evaluation (ggplot aes, parent frame vars)
utils::globalVariables(c(
  "model", "value", "Metric",  # BreakDiagnostic ggplot aes
  "loglike", "logmarglike", "sigma2start",  # BridgeChangeRegHybrid / BridgeMixedPanel
  "lambdadraws", "ns", "K", "sigmadraws", "alphadraws", "taudraws",  # marglike_BridgeChangeReg
  "sdraws", "n.break", "Pmat", "ntime", "ydm", "mu.st.state",
  "nu.shape", "nu.rate", "betadraws", "X", "c0", "d0",
  "known.alpha", "alpha.limit", "alpha.MH", "m", "A0", "y", "ps",
  "true.beta1", "true.beta2",  # dgp.break
  "gen_cov", "gen_beta_int", "rtnorm"  # sim / draw.beta (unused helper functions)
))


#' Fit Bayesian Bridge model with Change-points.
#'
#' This package provides Bayesian implementation of Bridge regression with possible change-points.
#' The package supports univariate time-series data as well as panel data with linear and non-linear link functions.
#' In addition to the Gaussian linear models, count and binary outcomes are supported.
#'
#' Specifically, this package offers the following functions:
#' \itemize{
#'   \item \code{\link{BridgeChangeReg}}: Univariate Regression with Change-points.
#'   \item \code{\link{BridgeFixedPanel}}: Fixed Effect Panel Regression with Change-points.
#' }
"_PACKAGE"
