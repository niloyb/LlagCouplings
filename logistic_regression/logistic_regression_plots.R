# devtools::install_github('cran/BayesLogit')
library(BayesLogit)
library(debiasedmcmc)
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

# Loading coupling times
source("tv_wasserstein_bounds.R")


# HMC small data times
hmc_meetings_lag_1000_hmc_stepsize_0.025_nsteps4 <- read.csv(file="hmc_meetings_lag_1000_stepsize_0.025_nsteps4.csv", header=TRUE, sep=",")
hmc_meetings_lag_1000_hmc_stepsize_0.025_nsteps5 <- read.csv(file="hmc_meetings_lag_1000_stepsize_0.025_nsteps5.csv", header=TRUE, sep=",")
hmc_meetings_lag_500_hmc_stepsize_0.025_nsteps6 <- read.csv(file="hmc_meetings_lag_500_stepsize_0.025_nsteps6.csv", header=TRUE, sep=",")
hmc_meetings_lag_2000_hmc_stepsize_0.025_nsteps7 <- read.csv(file="hmc_meetings_lag_2000_stepsize_0.025_nsteps7.csv", header=TRUE, sep=",")
# PG small data times
PG_meetings_lag_350 <- read.csv(file="PG_meetings_lag_350.csv", header=TRUE, sep=",")

t_start <- 0
t_end <- 2000
tv_ub_values_lag_1000_hmc_stepsize_0.025_nsteps4 <- sapply(c(t_start:t_end), function(t) tv_upper_bound_estimates(hmc_meetings_lag_1000_hmc_stepsize_0.025_nsteps4$x, L=1000, t))
tv_ub_values_lag_1000_hmc_stepsize_0.025_nsteps5 <- sapply(c(t_start:t_end), function(t) tv_upper_bound_estimates(hmc_meetings_lag_1000_hmc_stepsize_0.025_nsteps5$x, L=1000, t))
tv_ub_values_lag_500_hmc_stepsize_0.025_nsteps6 <- sapply(c(t_start:t_end), function(t) tv_upper_bound_estimates(hmc_meetings_lag_500_hmc_stepsize_0.025_nsteps6$x, L=500, t))
tv_ub_values_lag_2000_hmc_stepsize_0.025_nsteps7 <- sapply(c(t_start:t_end), function(t) tv_upper_bound_estimates(hmc_meetings_lag_2000_hmc_stepsize_0.025_nsteps7$x, L=2000, t))
tv_ub_values_PG_meetings_lag_350 <- sapply(c(t_start:t_end), function(t) tv_upper_bound_estimates(PG_meetings_lag_350$x, L=350, t))




######################## ggplot for bounds
ub_HMC_stepsize_stepsize_0.025_nsteps4 <- data.frame(cbind(c(t_start:t_end), colMeans(tv_ub_values_lag_1000_hmc_stepsize_0.025_nsteps4)))
colnames(ub_HMC_stepsize_stepsize_0.025_nsteps4) <- c('iterations', 'TV')
ub_HMC_stepsize_stepsize_0.025_nsteps5 <- data.frame(cbind(c(t_start:t_end), colMeans(tv_ub_values_lag_1000_hmc_stepsize_0.025_nsteps5)))
colnames(ub_HMC_stepsize_stepsize_0.025_nsteps5) <- c('iterations', 'TV')
ub_HMC_stepsize_stepsize_0.025_nsteps6 <- data.frame(cbind(c(t_start:t_end), colMeans(tv_ub_values_lag_500_hmc_stepsize_0.025_nsteps6)))
colnames(ub_HMC_stepsize_stepsize_0.025_nsteps6) <- c('iterations', 'TV')
ub_HMC_stepsize_stepsize_0.025_nsteps7 <- data.frame(cbind(c(t_start:t_end), colMeans(tv_ub_values_lag_2000_hmc_stepsize_0.025_nsteps7)))
colnames(ub_HMC_stepsize_stepsize_0.025_nsteps7) <- c('iterations', 'TV')
ub_PG <- data.frame(cbind(c(t_start:t_end), colMeans(tv_ub_values_PG_meetings_lag_350)))
colnames(ub_PG) <- c('iterations', 'TV')

g_tv_ub_HMC_PG <- ggplot(data=ub_PG,aes(x = iterations,y=TV, linetype = 'Polya-Gamma'))+geom_line(aes(y=TV)) + 
  geom_line(data=ub_HMC_stepsize_stepsize_0.025_nsteps4, aes(y=TV, linetype = 'L=4')) + 
  geom_line(data=ub_HMC_stepsize_stepsize_0.025_nsteps5, aes(y=TV, linetype = 'L=5')) +
  geom_line(data=ub_HMC_stepsize_stepsize_0.025_nsteps6, aes(y=TV, linetype = 'L=6')) +
  geom_line(data=ub_HMC_stepsize_stepsize_0.025_nsteps7, aes(y=TV, linetype = 'L=7')) +
  scale_x_continuous(expand = c(0, 0), limits = c(0,1800)) + ylim(0,1) +
  theme_grey(base_size = 18) + 
  theme(legend.position="bottom", legend.box = "horizontal") + labs(linetype="") + 
  scale_linetype_manual(values = c(c(2:5),1),
                        breaks=c("Polya-Gamma","L=4","L=5", "L=6","L=7"),
                        labels=c("Polya-Gamma", unname(TeX('$L_{HMC} = 0.1$')), unname(TeX('$L_{HMC} = 0.125$')),unname(TeX('$L_{HMC} = 0.15$')),unname(TeX('$L_{HMC} = 0.175$')))) + 
  labs(y = TeX("d_{TV}"))
g_tv_ub_HMC_PG 

ggsave(filename = "pg_vs_hmc_small_p.pdf", plot = g_tv_ub_HMC_PG, width = 9, height = 4.5)

# 



