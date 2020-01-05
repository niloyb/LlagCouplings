# This script reproduces Figure 5 in the article 
# "Estimating Convergence of Markov chains with L-Lag Couplings", 
# by Niloy Biswas, Pierre E. Jacob and Paul Vanetti. 
# In Advances in Neural Information Processing Systems 32 (2019), 7389--7399.
# https://papers.nips.cc/paper/8958-estimating-convergence-of-markov-chains-with-l-lag-couplings

# The data for the plots are generated using the ula_versus_mala.R script

#library(cowplot)
#library(dplyr)
#library(ggthemes)
#library(ggridges)
#library(grid)
#library(gridExtra)
#library(tidyverse)
#library(reshape2)
library(ggplot2)
library(latex2exp)

source("tv_wasserstein_bounds.R")


delta <- 0.25

df_list <- list()
for (dimension in c(50, 100, 200, 300, 400, 500, 600, 800, 1000)) {
  for (alg in c("ula", "mala")) {
    iteration_walltime <- read.csv(file=sprintf("%s_time_d=%d.csv", alg, dimension), header=TRUE, sep=",")$x
    
    meeting_times <- read.csv(file=sprintf("%s_meetings_d=%d_lag=10000.csv", alg, dimension), header=TRUE, sep=",")$V1

    t_start <- 0
    t_end <- 10000
    tv_ub_values <- sapply(t_start:t_end, function(t) tv_upper_bound_estimates(meeting_times, L=10000, t))

    mixing_time <- min(which(tv_ub_values < delta))
    adjusted_mixing_time <- mixing_time * iteration_walltime
    
    df_list[[length(df_list) + 1]] <- data.frame(
      alg=alg,
      dimension = dimension,
      mixing_time = adjusted_mixing_time
    )
    
  }
}
df <- bind_rows(df_list)

df["mixing_time"] <- df["mixing_time"] / min(df["mixing_time"])

g <- ggplot(data = df, aes(x = dimension, y = mixing_time, linetype = alg)) +
  geom_line() +
  xlab("dimension") + 
  ylab(TeX("$t_{mix}(0.25)$")) + 
  scale_x_continuous(expand = c(0, 0), limits = c(0,1020)) + 
  theme_grey(base_size = 16) +
  theme(legend.position = "bottom", legend.box = "horizontal") +
  labs(linetype = "")

ggsave(filename = "mala_vs_ula.pdf", plot = g, width = 8, height = 4)


