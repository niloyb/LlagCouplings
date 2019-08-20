# Bimodal target plots

#library(BayesLogit)
#library(debiasedmcmc)
library(latex2exp)
library(ggplot2)
library(grid)
library(gridExtra)
library(cowplot)
library(ggthemes)
library(ggridges)
library(reshape2)
library(tictoc)

library(dplyr)
# setmytheme()
rm(list = ls())
set.seed(21)
library(doParallel)
library(doRNG)
registerDoParallel(cores = detectCores()-2)

# Loading coupling functions
source("tv_wasserstein_bounds.R")

################################# Histogram of the 500th marginal
bad_bimodal_all_chains <- read.csv(file="bad_bimodal_traceplots.csv", header=TRUE, sep=",")
t_500_marginal <- bad_bimodal_all_chains[500,c(2:1001)]
trace_plot_df <- data.frame(t(t_500_marginal))
colnames(trace_plot_df) <- c('X')
hist <- ggplot(trace_plot_df) + geom_histogram(aes(x = X, y = ..density..),binwidth = 0.25, fill = "white", color = "black")
hist

## target distribution
target <- function(x){
  evals <- log(0.5) + dnorm(x, mean = c(-4, 4), sd = 1, log = TRUE)
  return(max(evals) + log(sum(exp(evals - max(evals)))))
}
curve(sapply(x, function(v) exp(target(v))), from = -10, to = 10)

bimodal_pdf_df <- data.frame(cbind(seq(-10,10,length.out = 10000), sapply(seq(-10,10,length.out = 10000), function(v) exp(target(v)))))
colnames(bimodal_pdf_df) <- c('x', 'pdf')

hist <- hist + geom_line(data = bimodal_pdf_df, aes(x = x, y = pdf), linetype = 1)
hist <- hist + theme_grey(base_size = 18)
hist

################################################################################################ TV bounds
bad_bimodal_lag_1 <- read.csv(file="bad_bimodal_lag_1.csv", header=TRUE, sep=",")
bad_bimodal_lag_18000 <- read.csv(file="bad_bimodal_lag_18000.csv", header=TRUE, sep=",")

t_start <- 0
t_end <- 10000
tv_ub_values_lag_1_bad_bimodal <- sapply(c(t_start:t_end), function(t) tv_upper_bound_estimates(bad_bimodal_lag_1$x, L=1, t))
tv_ub_values_lag_18000_bad_bimodal <- sapply(c(t_start:t_end), function(t) tv_upper_bound_estimates(bad_bimodal_lag_18000$x, L=18000, t))

######################## ggplot for TV bounds
lag_1_ub_small <- data.frame(cbind(c(t_start:t_end), colMeans(tv_ub_values_lag_1_bad_bimodal[c(1:1000),])))
lag_1_ub_med_1 <- data.frame(cbind(c(t_start:t_end), colMeans(tv_ub_values_lag_1_bad_bimodal[c(1:5000),])))
lag_1_ub_full <- data.frame(cbind(c(t_start:t_end), colMeans(tv_ub_values_lag_1_bad_bimodal[c(1:10000),])))
colnames(lag_1_ub_small) <- colnames(lag_1_ub_med_1)  <- colnames(lag_1_ub_full) <- c('iteration', 'TV')
lag_18000_ub_small <- data.frame(cbind(c(t_start:t_end), colMeans(tv_ub_values_lag_18000_bad_bimodal[c(1:1000),])))
lag_18000_ub_med_1 <- data.frame(cbind(c(t_start:t_end), colMeans(tv_ub_values_lag_18000_bad_bimodal[c(1:5000),])))
lag_18000_ub_full <- data.frame(cbind(c(t_start:t_end), colMeans(tv_ub_values_lag_18000_bad_bimodal[c(1:10000),])))
colnames(lag_18000_ub_small) <- colnames(lag_18000_ub_med_1) <- colnames(lag_18000_ub_full) <- c('iteration', 'TV')

g_tv_ub_hard <- ggplot(data=lag_1_ub_small,aes(x = iteration,y=TV, linetype = 'lag_1_small'))+geom_line(aes(y=TV)) + 
  geom_line(data=lag_1_ub_med_1, aes(y=TV, linetype = 'lag_1_med_1')) + 
  geom_line(data=lag_1_ub_full, aes(y=TV, linetype = 'lag_1_full')) + 
  geom_line(data=lag_18000_ub_small, aes(y=TV, linetype = 'lag_18000_small')) +
  geom_line(data=lag_18000_ub_med_1, aes(y=TV, linetype = 'lag_18000_med_1')) + 
  geom_line(data=lag_18000_ub_full, aes(y=TV, linetype = 'lag_18000_full')) + 
  theme_grey(base_size = 18) + 
  scale_x_continuous(trans = "log10",breaks = c(10^c(1:4))) + ylim(0,2) + labs(linetype = "") + labs(x = "iteration") +
  theme(legend.position="right") +
  scale_linetype_manual(values = c(c(4:2),rep(1,3)),breaks=c("lag_1_small", "lag_1_med_1", "lag_1_full", "lag_18000_small","lag_18000_med_1","lag_18000_full"), 
                        labels=c(unname(TeX("$\\tau^{(1)}_{1:1000}$")), unname(TeX("$\\tau^{(1)}_{1:5000}$")),
                                 unname(TeX("$\\tau^{(1)}_{1:10000}$")), unname(TeX("$\\tau^{(18000)}_{1:1000}$")),
                                 unname(TeX("$\\tau^{(18000)}_{1:5000}$")), unname(TeX("$\\tau^{(18000)}_{1:10000}$")))) +
  labs(y = TeX("d_{TV}")) 
g_tv_ub_hard


################################################################################################ Combined plot
combined_plot <- plot_grid(hist, g_tv_ub_hard, nrow=1, axis = 'b', rel_widths = c(0.42, 0.58))
ggsave(filename = "bimodal_hard_3.pdf", plot = combined_plot, width = 8, height = 4)







