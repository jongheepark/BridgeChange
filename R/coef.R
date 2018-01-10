


#' Obtain posterior mean of coefficients
#'
#' @param out
#' \code{BridgeChange} objects. Typically output from \code{BridgeChangeReg}.
#' @export
#'
coef.BridgeChange <- function(out) {
  estB  <- summary(out)$statistics[,1]
  estB <- estB[grep("beta", names(estB))]
  return(estB)
}

summarize_coef <- function(out, X, Y) {
  estB  <- summary(out)$statistics[,1]
  estB <- estB[grep("beta", names(estB))]
  estB <- estB * sd(Y) / apply(X, 2, sd)
  return(estB)
}



# plot.coef <- function(out, scale.back = FALSE, pooled.beta = NULL, true.beta= NULL, varnames=NULL){
#     library(tidyr)
#     library(dplyr)
#     library(scales)
#     library(tibble)
#     m = attr(out, "m");
#     ns <- m + 1
#     y=attr(out, "y");
#     X=attr(out, "X")
#     y.sd=attr(out, "y.sd");
#     X.sd=attr(out, "X.sd")
#     K <- ncol(X)
#     b.out <- out[, 1:(ns*K)]
#     if(scale.back){
#         b.out <- b.out*y.sd/rep(X.sd, ns)
#     }
#     b.mean <- matrix(apply(b.out, 2, mean), m+1, K, byrow=TRUE)
#     b.low <- matrix(apply(b.out, 2, quantile, probs=0.025), m+1, K, byrow=TRUE)
#     b.high <- matrix(apply(b.out, 2, quantile, probs=0.975), m+1, K, byrow=TRUE)
#     b.dist <- apply(b.mean, 2, dist)
#     b.avg = apply(b.mean, 2, mean)
#     if(!is.null(varnames)){
#         b.names <- varnames
#     }else{
#         b.names <- paste0("beta", 1:K)
#     }
#     ##
#     if(!is.null(true.beta)){
#         true.b.dist <- apply(true.beta, 2, dist)
#         true.b.avg = apply(true.beta, 2, mean)
#
#         coef0 <- tibble(Parameter = as.factor(rep(b.names, each=ns+2)),
#                         Regime = c(as.character(rep(c(paste0("Regime", 1:ns), paste0("True Regime", 1:ns)), K))),
#                         value = c(rbind(b.mean, true.beta)), dist = c(sapply(1:K, function(i){c(rep(b.dist[i], 2), rep(true.b.dist[i], 2))})),
#                         low = c(rbind(b.low, matrix(0, 2, K))), high = c(rbind(b.high, matrix(0, 2, K))),
#                         avg = c(sapply(1:K, function(i){c(rep(b.avg[i], 2), rep(true.b.avg[i], 2))})))
#
#      ggplot2::ggplot(coef0, aes(x=value, y=reorder(Parameter, dist))) +
#             ##  scale_color_discrete(labels = c("Pre-break", "Post-break")) +
#             geom_line(aes(group = Parameter), alpha = .4) +
#             geom_point(aes(color = Regime), size = 4.5, alpha = .4) +
#             scale_color_manual(name="Type",
#                                labels = c("Regime 1", "Regime 2", "True Regime 1", "True Regime 2"),
#                                values = c("blue", "navy", "tomato1", "tomato4")) +
#             ## geom_line(data = coef.regime, aes(group = Parameter)) +
#             ## geom_point(data = coef0, aes(x=true, y=reorder(Parameter, dist), color = Regime), size = 6, alpha=.2) +
#             ## geom_point(data = true.coef.regime, aes(color = regime), alpha = .3, size = 6) +
#             geom_vline(xintercept = 0, color = "black", alpha = .6) +
#             labs(title = "Parameter Changes",
#                  subtitle = " ") +
#             theme_minimal() +
#             theme(axis.title = element_blank(),
#                   panel.grid.major.x = element_blank(),
#                   panel.grid.minor = element_blank(),
#                   legend.title = element_blank(),
#                   legend.justification = c(0, 1),
#                   legend.position = c(.21, 1.075),
#                   legend.background = element_blank(),
#                   legend.direction="horizontal",
#                   text = element_text(family = "sans"),
#                   plot.title = element_text(size = 20, margin = margin(b = 10)),
#                   plot.subtitle = element_text(size = 10, color = "darkslategrey", margin = margin(b = 25)),
#                   plot.caption = element_text(size = 8, margin = margin(t = 10), color = "grey70", hjust = 0))
#
#     }else if(!is.null(pooled.beta)){
#
#        coef0 <- tibble(Parameter = as.factor(rep(b.names, each=ns+1)),
#                        Regime = c(as.character(rep(c(paste0("Regime", 1:ns), "Pooled"), K))),
#                        value = c(rbind(b.mean, pooled.beta)), dist = c(sapply(1:K, function(i){c(rep(b.dist[i], 2), 0)})),
#                        low = c(rbind(b.low), rep(0, K)), high = c(rbind(b.high, rep(0, K))),
#                        avg = c(sapply(1:K, function(i){c(rep(b.avg[i], 2), 0)})))
#
#         ggplot2::ggplot(coef0, aes(x=value, y=reorder(Parameter, dist))) +
#             ##  scale_color_discrete(labels = c("Pre-break", "Post-break")) +
#             geom_line(aes(group = Parameter), alpha = .4) +
#             geom_point(aes(color = Regime), size = 4.5, alpha = .4) +
#             scale_color_manual(name="Type",
#                                labels = c("Pooled", "Regime 1", "Regime 2"),
#                                values = c("tomato3", "blue", "navy")) +
#             ## geom_line(data = coef.regime, aes(group = Parameter)) +
#             ## geom_point(data = coef0, aes(x=true, y=reorder(Parameter, dist), color = Regime), size = 6, alpha=.2) +
#             ## geom_point(data = true.coef.regime, aes(color = regime), alpha = .3, size = 6) +
#             geom_vline(xintercept = 0, color = "black", alpha = .6) +
#             labs(title = "Parameter Changes",
#                  subtitle = " ") +
#             theme_minimal() +
#             theme(axis.title = element_blank(),
#                   panel.grid.major.x = element_blank(),
#                   panel.grid.minor = element_blank(),
#                   legend.title = element_blank(),
#                   legend.justification = c(0, 1),
#                   legend.position = c(.21, 1.075),
#                   legend.background = element_blank(),
#                   legend.direction="horizontal",
#                   text = element_text(family = "sans"),
#                   plot.title = element_text(size = 20, margin = margin(b = 10)),
#                   plot.subtitle = element_text(size = 10, color = "darkslategrey", margin = margin(b = 25)),
#                   plot.caption = element_text(size = 8, margin = margin(t = 10), color = "grey70", hjust = 0))
#
#
#     }else{
#         coef0 <- tibble(Parameter = as.factor(rep(b.names, each=ns)),
#                         regime = as.character(rep(paste0("Regime", 1:ns), K)),
#                         value = c(b.mean), dist = rep(c(b.dist), each =ns),
#                         low = c(b.low), high = c(b.high),
#                         avg = rep(c(b.avg), each =ns))
#
#         ggplot2::ggplot(coef0, aes(x=value, y=reorder(Parameter, dist))) +
#             ##  scale_color_discrete(labels = c("Pre-break", "Post-break")) +
#             geom_line(aes(group = Parameter), alpha = .8) +
#                 geom_point(aes(color = regime), size = 1.5, alpha = .8) +
#                     ## geom_line(data = coef.regime, aes(group = Parameter)) +
#                     ## geom_point(data = coef.regime, aes(color = regime), size = 2) +
#                     ## geom_point(data = true.coef.regime, aes(color = regime), alpha = .3, size = 6) +
#                     geom_vline(xintercept = 0, color = "black", alpha = .6) +
#                         labs(title = "Parameter Changes",
#                              subtitle = " ") +
#                                  theme_minimal() +
#                                      theme(axis.title = element_blank(),
#                                            panel.grid.major.x = element_blank(),
#                                            panel.grid.minor = element_blank(),
#                                            legend.title = element_blank(),
#                                            legend.justification = c(0, 1),
#                                            legend.position = c(.21, 1.075),
#                                            legend.background = element_blank(),
#                                            legend.direction="horizontal",
#                                            text = element_text(family = "sans"),
#                                            plot.title = element_text(size = 20, margin = margin(b = 10)),
#                                            plot.subtitle = element_text(size = 10, color = "darkslategrey", margin = margin(b = 25)),
#                                            plot.caption = element_text(size = 8, margin = margin(t = 10), color = "grey70", hjust = 0))
#
#     }
# }
