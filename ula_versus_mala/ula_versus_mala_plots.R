# ULA v MALA plots

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

source("tv_wasserstein_bounds.R")



################################################################################################ Dim 50
# ULA/ MALA meetings
ula_meetings_dim50_lag_1000 <- read.csv(file="ula_meetings_dim50_lag_1000.csv", header=TRUE, sep=",")
mala_meetings_dim50_lag_1000 <- read.csv(file="mala_meetings_dim50_lag_1000.csv", header=TRUE, sep=",")


t_start <- 0
t_end <- 10000
tv_ub_values_lag_1000_ula_dim50 <- sapply(c(t_start:t_end), function(t) tv_upper_bound_estimates(ula_meetings_dim50_lag_1000$x, L=1000, t))
tv_ub_values_lag_1000_mala_dim50 <- sapply(c(t_start:t_end), function(t) tv_upper_bound_estimates(mala_meetings_dim50_lag_1000$x, L=1000, t))


######################## ggplot for bounds
ub_ula_dim50 <- data.frame(cbind(c(t_start:t_end), colMeans(tv_ub_values_lag_1000_ula_dim50)))
ub_mala_dim50 <- data.frame(cbind(c(t_start:t_end), colMeans(tv_ub_values_lag_1000_mala_dim50)))
colnames(ub_ula_dim50) <- colnames(ub_mala_dim50) <- c('iterations', 'TV')

elapsedtime_ula_50 <- 0.000106
elapsedtime_mala_50 <- 0.000215


################################################################################################ Dim 100
# ULA/ MALA meetings
ula_meetings_dim100_lag_1000 <- read.csv(file="ula_meetings_dim100_lag_1000.csv", header=TRUE, sep=",")
mala_meetings_dim100_lag_1000 <- read.csv(file="mala_meetings_dim100_lag_1000.csv", header=TRUE, sep=",")

t_start <- 0
t_end <- 10000
tv_ub_values_lag_1000_ula_dim100 <- sapply(c(t_start:t_end), function(t) tv_upper_bound_estimates(ula_meetings_dim100_lag_1000$x, L=1000, t))
tv_ub_values_lag_1000_mala_dim100 <- sapply(c(t_start:t_end), function(t) tv_upper_bound_estimates(mala_meetings_dim100_lag_1000$x, L=1000, t))

######################## ggplot for bounds
ub_ula_dim100 <- data.frame(cbind(c(t_start:t_end), colMeans(tv_ub_values_lag_1000_ula_dim100)))
ub_mala_dim100 <- data.frame(cbind(c(t_start:t_end), colMeans(tv_ub_values_lag_1000_mala_dim100)))
colnames(ub_ula_dim100) <- colnames(ub_mala_dim100) <- c('iterations', 'TV')

elapsedtime_ula_100 <- 0.00025943
elapsedtime_mala_100 <- 0.0004863



################################################################################################ Dim 200
# ULA/ MALA meetings
ula_meetings_dim200_lag_1000 <- read.csv(file="ula_meetings_dim200_lag_1000.csv", header=TRUE, sep=",")
mala_meetings_dim200_lag_1000 <- read.csv(file="mala_meetings_dim200_lag_1000.csv", header=TRUE, sep=",")

t_start <- 0
t_end <- 10000
tv_ub_values_lag_1000_ula_dim200 <- sapply(c(t_start:t_end), function(t) tv_upper_bound_estimates(ula_meetings_dim200_lag_1000$x, L=1000, t))
tv_ub_values_lag_1000_mala_dim200 <- sapply(c(t_start:t_end), function(t) tv_upper_bound_estimates(mala_meetings_dim200_lag_1000$x, L=1000, t))

######################## ggplot for bounds
ub_ula_dim200 <- data.frame(cbind(c(t_start:t_end), colMeans(tv_ub_values_lag_1000_ula_dim200)))
ub_mala_dim200 <- data.frame(cbind(c(t_start:t_end), colMeans(tv_ub_values_lag_1000_mala_dim200)))
colnames(ub_ula_dim200) <- colnames(ub_mala_dim200) <- c('iterations', 'TV')

elapsedtime_ula_200 <-  0.0008
elapsedtime_mala_200 <- 0.001818

g_tv_ub_ula_mala_200 <- ggplot(data=ub_ula_dim200,aes(x = iterations,y=TV, linetype = 'ULA'))+geom_line(aes(y=TV)) + 
  geom_line(data=ub_mala_dim200, aes(x = iterations, y=TV, linetype = 'MALA')) +
  geom_line(data=ub_ula_dim200, aes(x = iterations*(elapsedtime_mala_200/elapsedtime_ula_200), y=TV, linetype = 'Time-adjusted MALA')) +
  scale_x_continuous(expand = c(0, 0), limits = c(0,2200)) + ylim(0,1) +
  theme_grey(base_size = 20) + 
  theme(legend.position="bottom", legend.box = "horizontal") + labs(linetype="Chain") + 
  labs(y = TeX("d_{TV}"))
g_tv_ub_ula_mala_200 
#

################################################################################################ Dim 300
# ULA/ MALA meetings
ula_meetings_dim300_lag_2000 <- read.csv(file="ula_meetings_dim300_lag_2000.csv", header=TRUE, sep=",")
mala_meetings_dim300_lag_2000 <- read.csv(file="mala_meetings_dim300_lag_2000.csv", header=TRUE, sep=",")

t_start <- 0
t_end <- 10000
tv_ub_values_lag_2000_ula_dim300 <- sapply(c(t_start:t_end), function(t) tv_upper_bound_estimates(ula_meetings_dim300_lag_2000$x, L=2000, t))
tv_ub_values_lag_2000_mala_dim300 <- sapply(c(t_start:t_end), function(t) tv_upper_bound_estimates(mala_meetings_dim300_lag_2000$x, L=2000, t))

######################## ggplot for bounds
ub_ula_dim300 <- data.frame(cbind(c(t_start:t_end), colMeans(tv_ub_values_lag_2000_ula_dim300)))
ub_mala_dim300 <- data.frame(cbind(c(t_start:t_end), colMeans(tv_ub_values_lag_2000_mala_dim300)))
colnames(ub_ula_dim300) <- colnames(ub_mala_dim300) <- c('iterations', 'TV')

elapsedtime_ula_300 <- 0.001947
elapsedtime_mala_300 <- 0.003619


################################################################################################ Dim 400
# RWMH/ ULA/ MALA meetings
ula_meetings_dim400_lag_2000 <- read.csv(file="ula_meetings_dim400_lag_2000.csv", header=TRUE, sep=",")
mala_meetings_dim400_lag_2000 <- read.csv(file="mala_meetings_dim400_lag_2000.csv", header=TRUE, sep=",")

t_start <- 0
t_end <- 10000
tv_ub_values_lag_2000_ula_dim400 <- sapply(c(t_start:t_end), function(t) tv_upper_bound_estimates(ula_meetings_dim400_lag_2000$x, L=2000, t))
tv_ub_values_lag_2000_mala_dim400 <- sapply(c(t_start:t_end), function(t) tv_upper_bound_estimates(mala_meetings_dim400_lag_2000$x, L=2000, t))

######################## ggplot for bounds
ub_ula_dim400 <- data.frame(cbind(c(t_start:t_end), colMeans(tv_ub_values_lag_2000_ula_dim400)))
ub_mala_dim400 <- data.frame(cbind(c(t_start:t_end), colMeans(tv_ub_values_lag_2000_mala_dim400)))
colnames(ub_ula_dim400) <- colnames(ub_mala_dim400) <- c('iterations', 'TV')

elapsedtime_ula_400 <-  0.00262
elapsedtime_mala_400 <- 0.006003

################################################################################################ Dim 500
# RWMH/ ULA/ MALA meetings
ula_meetings_dim500_lag_2000 <- read.csv(file="ula_meetings_dim500_lag_2000.csv", header=TRUE, sep=",")
mala_meetings_dim500_lag_4000 <- read.csv(file="mala_meetings_dim500_lag_4000.csv", header=TRUE, sep=",")
t_start <- 0
t_end <- 10000
tv_ub_values_lag_2000_ula_dim500 <- sapply(c(t_start:t_end), function(t) tv_upper_bound_estimates(ula_meetings_dim500_lag_2000$x, L=2000, t))
tv_ub_values_lag_4000_mala_dim500 <- sapply(c(t_start:t_end), function(t) tv_upper_bound_estimates(mala_meetings_dim500_lag_4000$x, L=4000, t))
######################## ggplot for bounds
ub_ula_dim500 <- data.frame(cbind(c(t_start:t_end), colMeans(tv_ub_values_lag_2000_ula_dim500)))
ub_mala_dim500 <- data.frame(cbind(c(t_start:t_end), colMeans(tv_ub_values_lag_4000_mala_dim500)))
colnames(ub_ula_dim500) <- colnames(ub_mala_dim500) <- c('iterations', 'TV')
elapsedtime_ula_500 <-  0.003798
elapsedtime_mala_500 <- 0.009442

################################################################################################ Dim 600
# RWMH/ ULA/ MALA meetings
ula_meetings_dim600_lag_2000 <- read.csv(file="ula_meetings_dim600_lag_2000.csv", header=TRUE, sep=",")
mala_meetings_dim600_lag_4000 <- read.csv(file="mala_meetings_dim600_lag_4000.csv", header=TRUE, sep=",")
t_start <- 0
t_end <- 10000
tv_ub_values_lag_2000_ula_dim600 <- sapply(c(t_start:t_end), function(t) tv_upper_bound_estimates(ula_meetings_dim600_lag_2000$x, L=2000, t))
tv_ub_values_lag_4000_mala_dim600 <- sapply(c(t_start:t_end), function(t) tv_upper_bound_estimates(mala_meetings_dim600_lag_4000$x, L=4000, t))
######################## ggplot for bounds
ub_ula_dim600 <- data.frame(cbind(c(t_start:t_end), colMeans(tv_ub_values_lag_2000_ula_dim600)))
ub_mala_dim600 <- data.frame(cbind(c(t_start:t_end), colMeans(tv_ub_values_lag_4000_mala_dim600)))
colnames(ub_ula_dim600) <- colnames(ub_mala_dim600) <- c('iterations', 'TV')
elapsedtime_ula_600 <-  0.008281
elapsedtime_mala_600 <- 0.016414

################################################################################################ Dim 800
# RWMH/ ULA/ MALA meetings
ula_meetings_dim800_lag_10000 <- read.csv(file="ula_meetings_dim800_lag_10000.csv", header=TRUE, sep=",")
mala_meetings_dim800_lag_10000 <- read.csv(file="mala_meetings_dim800_lag_10000.csv", header=TRUE, sep=",")

t_start <- 0
t_end <- 10000
tv_ub_values_lag_10000_ula_dim800 <- sapply(c(t_start:t_end), function(t) tv_upper_bound_estimates(ula_meetings_dim800_lag_10000$x, L=10000, t))
tv_ub_values_lag_10000_mala_dim800 <- sapply(c(t_start:t_end), function(t) tv_upper_bound_estimates(mala_meetings_dim800_lag_10000$x, L=10000, t))
######################## ggplot for bounds
ub_ula_dim800 <- data.frame(cbind(c(t_start:t_end), colMeans(tv_ub_values_lag_10000_ula_dim800)))
ub_mala_dim800 <- data.frame(cbind(c(t_start:t_end), colMeans(tv_ub_values_lag_10000_mala_dim800)))
colnames(ub_ula_dim800) <- colnames(ub_mala_dim800) <- c('iterations', 'TV')
elapsedtime_ula_800 <-  0.01083
elapsedtime_mala_800 <- 0.0253


################################################################################################ Dim 1000
# RWMH/ ULA/ MALA meetings
ula_meetings_dim1000_lag_10000 <- read.csv(file="ula_meetings_dim1000_lag_10000.csv", header=TRUE, sep=",")
mala_meetings_dim1000_lag_10000 <- read.csv(file="mala_meetings_dim1000_lag_10000.csv", header=TRUE, sep=",")
t_start <- 0
t_end <- 10000
tv_ub_values_lag_10000_ula_dim1000 <- sapply(c(t_start:t_end), function(t) tv_upper_bound_estimates(ula_meetings_dim1000_lag_10000$x, L=10000, t))
tv_ub_values_lag_10000_mala_dim1000 <- sapply(c(t_start:t_end), function(t) tv_upper_bound_estimates(mala_meetings_dim1000_lag_10000$x, L=10000, t))
######################## ggplot for bounds
ub_ula_dim1000 <- data.frame(cbind(c(t_start:t_end), colMeans(tv_ub_values_lag_10000_ula_dim1000)))
ub_mala_dim1000 <- data.frame(cbind(c(t_start:t_end), colMeans(tv_ub_values_lag_10000_mala_dim1000)))
colnames(ub_ula_dim1000) <- colnames(ub_mala_dim1000) <- c('iterations', 'TV')
elapsedtime_ula_1000 <-  0.01804
elapsedtime_mala_1000 <- 0.04517


################################################################################################ Mixing times plot
delta <- 0.25
mixing_time <- function(tv_bound_values, delta){ return(min(which(tv_bound_values<delta))) }
ula_mixing <- c(mixing_time(ub_ula_dim50$TV, delta), mixing_time(ub_ula_dim100$TV, delta), mixing_time(ub_ula_dim200$TV, delta), mixing_time(ub_ula_dim300$TV, delta), mixing_time(ub_ula_dim400$TV, delta), mixing_time(ub_ula_dim500$TV, delta), mixing_time(ub_ula_dim600$TV, delta), mixing_time(ub_ula_dim800$TV, delta), mixing_time(ub_ula_dim1000$TV, delta))
mala_mixing <- c(mixing_time(ub_mala_dim50$TV, delta), mixing_time(ub_mala_dim100$TV, delta), mixing_time(ub_mala_dim200$TV, delta), mixing_time(ub_mala_dim300$TV, delta), mixing_time(ub_mala_dim400$TV, delta), mixing_time(ub_mala_dim500$TV, delta), mixing_time(ub_mala_dim600$TV, delta), mixing_time(ub_mala_dim800$TV, delta), mixing_time(ub_mala_dim1000$TV, delta))
time_adjutment <- c(elapsedtime_mala_50/elapsedtime_ula_50,elapsedtime_mala_100/elapsedtime_ula_100,elapsedtime_mala_200/elapsedtime_ula_200,elapsedtime_mala_300/elapsedtime_ula_300,elapsedtime_mala_400/elapsedtime_ula_400,elapsedtime_mala_500/elapsedtime_ula_500, elapsedtime_mala_600/elapsedtime_ula_600, elapsedtime_mala_800/elapsedtime_ula_800, elapsedtime_mala_1000/elapsedtime_ula_1000)
mala_time_adjusted_mixing <- mala_mixing*time_adjutment
mixing_time_w_iteration_df <- data.frame(cbind(c(50,100,200,300,400,500,600,800,1000), ula_mixing, mala_mixing))
colnames(mixing_time_w_iteration_df) <- c('p', 'ULA', 'MALA')

mixing_time_w_iteration_df_melted <- melt(mixing_time_w_iteration_df, id.vars = 'p')
colnames(mixing_time_w_iteration_df_melted) <- c('p', 'chain_type', 't_mixing')


g_1 <- ggplot(data=mixing_time_w_iteration_df_melted,aes(x=p, y = t_mixing, linetype=chain_type)) + geom_line() 
g_1 <- g_1 + xlab("dimension") + ylab(TeX("$t_{mix}(0.25)$")) + 
  scale_x_continuous(expand = c(0, 0), limits = c(0,1020)) + theme_grey(base_size = 16) +
  theme(legend.position="bottom", legend.box = "horizontal") + labs(linetype="")
g_1

#ggsave(filename = "mala_vs_ula_2.pdf", plot = g_1, width = 8, height = 4)


