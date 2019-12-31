# This script reproduces Figure 4 in the article 
# "Estimating Convergence of Markov chains with L-Lag Couplings", 
# by Niloy Biswas, Pierre E. Jacob and Paul Vanetti. 
# In Advances in Neural Information Processing Systems 32 (2019), 7389--7399.
# https://papers.nips.cc/paper/8958-estimating-convergence-of-markov-chains-with-l-lag-couplings

# The data for the plots are generated using the 
# polya_gamma.R and hmc_logistic_regression.R scripts

library(latex2exp)
library(tidyverse)

source("tv_wasserstein_bounds.R")

hmc_experiments_df <- data.frame(
  stepsize = c(0.025, 0.025, 0.025, 0.025),
  nsteps = c(4, 5, 6, 7),
  lag = c(1000, 1000, 500, 2000)
)

t_start <- 0
t_end <- 2000

bound_df_list <- list()
for (i in 1:nrow(hmc_experiments_df)) {
  row <- hmc_experiments_df[i, ]
  filename = sprintf("hmc_meetings_lag=%d_stepsize=%f_nsteps=%d.csv", row$lag, row$stepsize, row$nsteps)
  meetings <- read.csv(file=filename, header=TRUE, sep=",")
  tv_bounds <- colMeans(sapply(c(t_start:t_end), function(t) tv_upper_bound_estimates(meetings$x, L=row$lag, t)))
  bound_df_list[[length(bound_df_list) + 1]] <- data.frame(t=t_start:t_end, bounds=tv_bounds, alg=sprintf("hmc_L=%d", row$nsteps))
}

# Polya-Gamma meeting times.
pg_meetings_lag_350 <- read.csv(file="PG_meetings_lag=350.csv", header=TRUE, sep=",")
tv_bounds <- colMeans(sapply(c(t_start:t_end), function(t) tv_upper_bound_estimates(pg_meetings_lag_350$x, L=350, t)))
bound_df_list[[length(bound_df_list) + 1]] <- data.frame(t=t_start:t_end, bounds=tv_bounds, alg="pg")

bound_df <- bind_rows(bound_df_list)

pg_vs_hmc_plot <- ggplot(data=bound_df, aes(x = t, y=bounds, linetype = alg)) +
  geom_line() + 
  scale_x_continuous(expand = c(0, 0), limits = c(0,1800)) + 
  scale_linetype_manual(
    values = c(c(2:10),1),
    breaks=c("pg", "hmc_L=4", "hmc_L=5", "hmc_L=6", "hmc_L=7"),
    labels=c(
      "Polya-Gamma",
      unname(TeX('$L_{HMC} = 0.1$')),
      unname(TeX('$L_{HMC} = 0.125$')),
      unname(TeX('$L_{HMC} = 0.15$')),
      unname(TeX('$L_{HMC} = 0.175$'))
    )) +
  ylim(0,1) +
  theme_grey(base_size = 18) + 
  theme(legend.position="bottom", legend.box = "horizontal") +
  labs(linetype="") + 
  labs(y = TeX("d_{TV}"))

ggsave(filename = "pg_vs_hmc_small_p.pdf", plot = pg_vs_hmc_plot, width = 9, height = 4.5)
