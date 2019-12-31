# This script reproduces Figure 2 in the article 
# "Estimating Convergence of Markov chains with L-Lag Couplings", 
# by Niloy Biswas, Pierre E. Jacob and Paul Vanetti. 
# In Advances in Neural Information Processing Systems 32 (2019), 7389--7399.
# https://papers.nips.cc/paper/8958-estimating-convergence-of-markov-chains-with-l-lag-couplings

# The data for the plots are generated using the bimodal_target.R script

library(cowplot)
library(dplyr)
library(ggplot2)
library(latex2exp)
library(purrr)

source("tv_wasserstein_bounds.R")

# Histogram of the 500th marginal.
bad_bimodal_500_df <- read.csv(file="bad_bimodal_500.csv", header=TRUE, sep=",")

# PDF of target distribution.
logtarget <- function(x){
  evals <- log(0.5) + dnorm(x, mean = c(-4, 4), sd = 1, log = TRUE)
  return(max(evals) + log(sum(exp(evals - max(evals)))))
}
bimodal_pdf_df <- data.frame(x = seq(-10, 10, length.out = 10000),
                             pdf = sapply(seq(-10, 10, length.out = 10000), function(v) exp(logtarget(v))))

# Compute TV bounds from meeting times.
t_start <- 0
t_end <- 10000

lags <- c(1, 18000)
ns <- c(1000, 5000, 10000)

df_list <- list()
for (lag in lags) {
  filename <- sprintf("bad_bimodal_lag=%d.csv", lag)
  meeting_times <- read.csv(file=filename, header=TRUE, sep=",")
  tv_ub_values <- sapply(t_start:t_end, function(t) tv_upper_bound_estimates(meeting_times$V1, L = lag, t))
  for (nsamples in ns) {
    df <- data.frame(t=c(t_start:t_end), tv_bound=colMeans(tv_ub_values[1:nsamples, ]), nsamples=nsamples, lag=lag)
    df_list[[length(df_list) + 1]] <- df
  }
}
df <- bind_rows(df_list)

hist <- ggplot(bad_bimodal_500_df) +
  geom_histogram(aes(x = V1, y = ..density..), binwidth = 0.25, fill = "white", color = "black") +
  geom_line(data = bimodal_pdf_df, aes(x = x, y = pdf), linetype = 1) +
  theme_grey(base_size = 18)
hist

g_tv_ub_hard <- ggplot(data = df, aes(x = t, y = tv_bound, linetype = paste(lag, nsamples))) +
  geom_line() +
  theme_grey(base_size = 18) + 
  scale_x_continuous(trans = "log10", breaks = c(10 ^ c(1:4))) +
  ylim(0,2) +
  labs(linetype = "") +
  labs(x = "iteration") +
  theme(legend.position="right") +
  scale_linetype_manual(values = c(c(4:2), rep(1,3)),
                        breaks = sapply(cross2(ns, lags), function(x) sprintf("%d %d", x[[2]], x[[1]])),
                        labels = sapply(cross2(ns, lags), function(x) unname(TeX(sprintf("$\\tau^{(%d)}_{1:%d}$", x[[2]], x[[1]]))))) +
  labs(y = TeX("d_{TV}")) 
g_tv_ub_hard

combined_plot <- plot_grid(hist, g_tv_ub_hard, nrow=1, axis = 'b', rel_widths = c(0.42, 0.58))
ggsave(filename = "bimodal.pdf", plot = combined_plot, width = 8, height = 4)







