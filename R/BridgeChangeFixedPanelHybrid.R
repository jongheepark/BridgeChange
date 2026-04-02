######################################################################################################
## The idea is to develop a sparsity-induced prior-posterior model that fits data
## Bridge regression using mixture of normals representation.
## by JHP "Wed Oct 19 15:13:47 2016"
######################################################################################################

#' Adaptive Lasso Variable Selection
#'
#' Performs adaptive lasso variable selection using ridge regression for initial weight estimation.
#'
#' @param y Response vector.
#' @param x Design matrix.
#'
#' @return A vector of adaptive lasso coefficients.
#'
#' @export
#'
#'
## Adaptive for each regime
adaptive.lasso <- function(y, x) {
  fit.ridge <- cv.glmnet(
    y = y, x = x, type.measure = "mse",
    alpha = 0, standardize = TRUE, family = "gaussian"
  )
  w3 <- 1 / abs(matrix(coef(fit.ridge, s = fit.ridge$lambda.min)[, 1][2:(ncol(x) + 1)]))^1
  w3[w3[, 1] == Inf] <- 999999999

  cv.adaptive <- cv.glmnet(
    x = x, y = y, family = "gaussian",
    alpha = 1, penalty.factor = w3
  )

  beta.adaptive <- coef(cv.adaptive)[-1]
  return(beta.adaptive)
}

BridgeFixedPanelHybrid <- function(formula, data, index, model, effect,
                                   standardize = TRUE, interaction = 1,
                                   n.break = 1, ols.weight = FALSE,
                                   sparse.only = FALSE, alpha.MH = FALSE,
                                   mcmc = 100, burn = 100, verbose = 100,
                                   thin = 1, c0 = 0.1, d0 = 0.1,
                                   r0 = 1, R0 = 1, nu.shape = 2.0,
                                   nu.rate = 2.0, alpha = 1,
                                   Waic = FALSE, marginal = FALSE) {

  call <- match.call()
  a <- NULL
  b <- NULL

  prep <- .bc_prepare_panel(
    formula = formula,
    data    = data,
    index   = index,
    model   = model,
    effect  = effect
  )

  X <- prep$X
  y <- prep$y
  plm.index <- prep$idx

  ns <- n.break + 1
  plmX <- X
  plmy <- y
  m <- n.break
  W <- matrix(0, nrow = length(y), ncol = 1)
  colnames(W) <- "fixed_placeholder"

  interaction <- min(ncol(plmX), interaction)

  if (interaction > 1) {
    newX <- list()
    newX[[1]] <- X
    x1.1 <- data.frame(X, stringsAsFactors = FALSE)
    var.names <- colnames(X)

    for (j in 2:interaction) {
      x1.2 <- matrix(
        t(apply(x1.1, 1, combn, j, prod)),
        nrow = nrow(X)
      )
      newX[[j]] <- as.matrix(x1.2)
      colnames(newX[[j]]) <- combn(var.names, j, paste, collapse = "-")
    }

    X <- Reduce(cbind, newX)

    if (sum(apply(X, 2, sd) == 0) > 0) {
      cat("Some interactions (", sum(apply(X, 2, sd) == 0), ") are all zero. So they are removed!\n")
      X <- X[, apply(X, 2, sd) != 0, drop = FALSE]
    }
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

  var.names <- colnames(X)
  subject.id <- as.numeric(as.factor(plm.index[[1]]))
  time.id    <- as.numeric(as.factor(plm.index[[2]]))

  if (sparse.only) {
    output <- NA
  } else {
    output <- BridgeMixedPanel(
      subject.id = subject.id,
      time.id    = time.id,
      y          = as.vector(y),
      X          = X,
      W          = W,
      n.break    = n.break,
      c0         = c0,
      d0         = d0,
      r0         = r0,
      R0         = R0,
      standardize = FALSE,
      alpha.MH   = alpha.MH,
      mcmc       = mcmc,
      burn       = burn,
      thin       = thin,
      verbose    = verbose,
      nu.shape   = nu.shape,
      nu.rate    = nu.rate,
      alpha      = alpha,
      Waic       = Waic,
      marginal   = marginal,
      fixed      = TRUE,
      unscaled.Y = unscaled.Y,
      unscaled.X = unscaled.X
    )
  }

  ninv.y <- 1 / length(y)

  if (n.break > 0) {
    state <- round(apply(attr(output, "s.store"), 2, mean))
    cat("estiamted states are ", table(state), "\n")

    beta.mean <- matrix(
      apply(output, 2, mean)[grepl("beta", colnames(output))],
      ncol(X), ns
    )

    s2.mean <- matrix(
      apply(output, 2, mean)[grepl("sigma", colnames(output))],
      , ns
    )

    yhat.state <- X %*% beta.mean
    state.indicator <- state[time.id]
  }

  attr(output, "title") <- "BridgeChangeFixedPanelHybrid Posterior Sample"
  attr(output, "m") <- n.break
  attr(output, "plm_index") <- plm.index

  if (standardize) {
    attr(output, "dat.sd") <- dat.sd
  }

  return(output)
}