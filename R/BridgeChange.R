#' @useDynLib BridgeChange
#' @import RcppArmadillo
#' @importFrom Rcpp sourceCpp
#' @exportPattern "^[[:alpha:]]+"
NULL
#> NULL


#' Fit Bayesian Bridge model with Change-points.
#'
#' This package provides Bayesian implementation of Bridge regression with possible change-points.
#' The package supports univariate time-series data as well as panel data with linear and non-linear link functions.
#' In addition to the Gaussian linear models, count and binary outcomes are supported.
#'
#' Specifically, this package offers the following functions:
#' \itemize{
#'   \item \code{\link{BridgeChangeReg}}: Univariate Regression with Change-points.
#'   \item \code{\link{BridgeMixedPanel}}: Mixed Effect Panel Regression with Change-points.
#'   \item \code{\link{BridgeFixedPanel}}: Fixed Effect Panel Regression with Change-points.
#' }
"_PACKAGE"
