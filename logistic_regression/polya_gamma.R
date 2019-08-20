
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

source("tv_wasserstein_bounds.R")

## Loading data
# load german credit dataset
data(germancredit)
X <- scale(X)
X[,1] <- rep(1, nrow(X))
n <- nrow(X)
p <- ncol(X)

# Prior
b <- matrix(0, nrow = p, ncol = 1)
B <- diag(10, p, p)
logistic_setting <- logisticregression_precomputation(Y, X, b, B)

# Functions
single_kernel <- function(chain_state, logistic_setting){
  zs <- abs(logisticregression_xbeta(logistic_setting$X, t(chain_state)))
  w <- BayesLogit::rpg(logistic_setting$n, h=1, z=zs)
  res <- logisticregression_m_and_sigma(w, X, logistic_setting$invB, logistic_setting$KTkappaplusinvBtimesb)
  chain_state <- t(fast_rmvnorm_chol(1, res$m, res$Cholesky))
  return(chain_state)
}

coupled_kernel <- function(chain_state1, chain_state2, logistic_setting, return_ws=FALSE){
  ws <- sample_w(chain_state1, chain_state2, logistic_setting$X)
  betas <- sample_beta(ws$w1, ws$w2, logistic_setting)
  if (all(ws$w1 == ws$w2)){
    betas$beta2 <- betas$beta1
  }
  if(!return_ws){
    return(list(chain_state1=cbind(betas$beta1), chain_state2=cbind(betas$beta2)))
  } else {
    return(list(chain_state1=cbind(betas$beta1), chain_state2=cbind(betas$beta2),w1=ws$w1,w2=ws$w2))
  }
}

rinit <- function(){
  t(fast_rmvnorm(1, mean = b, covariance = B))
}

# Behavior from a single run
# current_value1 <- rinit()
# current_value2 <- rinit()
# niterations <- 1000
# 
# chain1_w <- chain2_w <- matrix(ncol=n, nrow=niterations)
# chain1_beta <- chain2_beta <- matrix(ncol=p, nrow=niterations)
# chain1_beta[1,] <- current_value1
# chain2_beta[1,] <- current_value2
# 
# ptc <- proc.time()
# for(t in 2:niterations){
#   current_value <- coupled_kernel(current_value1, current_value2, logistic_setting, return_ws=TRUE)
#   chain1_w[t-1,] <- current_value$w1
#   chain2_w[t-1,] <- current_value$w2
#   chain1_beta[t,] <- current_value1 <- current_value$chain_state1
#   chain2_beta[t,] <- current_value2 <- current_value$chain_state2
#   #if(all(current_value1==current_value2)) break
# }
# meetingtime=t
# time_taken <- (proc.time()[3]-ptc[3])/t
# 
# # Traceplot
# matplot(chain1_beta[c(100:1000),c(1:10)], type = 'l')
# 
# # gg Traceplot
# trace_plot_df <- data.frame(cbind(c(1:200),chain1_beta[c(1:200),c(1:49)]))
# trace_plot_df <- melt(trace_plot_df, id="X1")
# colnames(trace_plot_df) <- c('t', 'chain', 'X')
# g_PG <- ggplot(data=trace_plot_df,aes(x=t, y=X, colour=chain)) + geom_line() + theme_gray(base_size = 20)
# mcmc_traceplot_PG <- g_PG + scale_colour_grey(end = 0.5) + theme(legend.position = "none")
# mcmc_traceplot_PG

# Meeting times
polya_gamma_meeting_time <- function(single_kernel, coupled_kernel, rinit, max_iterations = Inf, L=1){
  
  if(max_iterations != Inf)
  {
    chain1 <- matrix(0, nrow = max(max_iterations,L), ncol = p)
    chain2 <- matrix(0, nrow = max(max_iterations,L), ncol = p)
  } else {
    chain1 <- matrix(0, nrow = 1000, ncol = p)
    chain2 <- matrix(0, nrow = 1000, ncol = p)
  }
  
  chain_state1 <- rinit()
  chain_state2 <- rinit()
  if (L>1)
  {
    for (l in 1:(L-1))
    {
      chain_state1 <- single_kernel(chain_state1, logistic_setting)
    }
  }
  
  chain1[1,] <- chain_state1
  chain2[1,] <- chain_state2
  
  
  chain1[2,] <- single_kernel(chain1[1,], logistic_setting)
  chain_state1 <- chain1[2,]
  
  iter <- L
  meet <- FALSE
  finished <- FALSE
  meetingtime <- Inf
  # iter here is L; at this point we have X_L,Y_0 and we are going to generate successively X_t,Y_{t-1} where iter = t
  while (!finished && iter < max_iterations){
    iter <- iter + 1
    if (meet){
      chain_state1 <- single_kernel(chain_state1)
      chain_state2 <- chain_state1
    } else {
      res_coupled_kernel <- coupled_kernel(chain_state1, chain_state2, logistic_setting)
      chain_state1 <- res_coupled_kernel$chain_state1
      chain_state2 <- res_coupled_kernel$chain_state2
      if (all(chain_state1 == chain_state2) && !meet){
        # recording meeting time tau
        meet <- TRUE
        meetingtime <- iter
      }
    }
    
    chain1[((iter)+2-L),] <- chain_state1
    chain2[((iter)+1-L),] <- chain_state2
    
    # Appending extra memory for chains
    if ((iter)+2-L+1000> dim(chain1)[1] )
    {
      chain1 <- rbind(chain1, matrix(0, nrow = 1000, ncol = p))
      chain2 <- rbind(chain2, matrix(0, nrow = 1000, ncol = p))
    }
    
    # stop after max(m, tau) steps
    if (iter >= meetingtime){
      finished <- TRUE
    }
  }
  return(list(meetingtime = meetingtime, chain1 = chain1, chain2 = chain2, iteration = iter, finished = finished))
}


# Generating meeting times
nsamples <- 1000
max_iterations <- 20000
t_start <- 0
t_end <- 1000
meetings_lag_350 <- rep(0, nsamples)
tv_ub_values_lag_350 <- matrix(0, nrow = nsamples, ncol = (t_end-t_start+1))
wass_ub_values_lag_350 <- matrix(0, nrow = nsamples, ncol = (t_end-t_start+1))

for(irep in 1:nsamples){
  lag_350 <- polya_gamma_meeting_time(single_kernel, coupled_kernel, rinit, max_iterations, L = 350)
  meetings_lag_350[irep] <- lag_350$meetingtime
  #tv_ub_values_lag_350[irep,] <- sapply(c(t_start:t_end), function(t) tv_upper_bound_estimates(lag_350$meetingtime, L=350, t))
  #wass_ub_values_lag_350[irep,] <- sapply(c(t_start:t_end), function(t) wasserstein_distance_upper_bound_estimate(lag_350$meetingtime, lag_350$chain1, lag_350$chain2, L=350, t))
  
  print(irep)
}

#write.csv(meetings_lag_350, "PG_meetings_lag_350.csv")
































upper_bound_estimates <- function(coupling_times, L, t)
{
  return(pmax(0,ceiling((coupling_times-L-t)/L)))
}

######################## ggplot for TV bounds
lag_1_tv_ub <- data.frame(cbind(c(1:200),sapply(1:200, function(t) mean(upper_bound_estimates(meetingtimes_1, 1, t))),
                                sapply(1:200, function(t) sd(upper_bound_estimates(meetingtimes_1, 1, t))/sqrt(nsamples))))
colnames(lag_1_tv_ub) <- c('t', 'means', 'sd')
lag_10_tv_ub <- data.frame(cbind(c(1:200),sapply(1:200, function(t) mean(upper_bound_estimates(meetingtimes_10, 10, t))),
                                 sapply(1:200, function(t) sd(upper_bound_estimates(meetingtimes_10, 10, t))/sqrt(nsamples))))
colnames(lag_10_tv_ub) <- c('t', 'means', 'sd')
lag_50_tv_ub <- data.frame(cbind(c(1:200),sapply(1:200, function(t) mean(upper_bound_estimates(meetingtimes_50, 50, t))),
                                 sapply(1:200, function(i) sd(upper_bound_estimates(meetingtimes_50, 50, t))/sqrt(nsamples))))
colnames(lag_50_tv_ub) <- c('t', 'means', 'sd')
lag_100_tv_ub <- data.frame(cbind(c(1:200),sapply(1:200, function(t) mean(upper_bound_estimates(meetingtimes_100, 100, t))),
                                  sapply(1:200, function(i) sd(upper_bound_estimates(meetingtimes_100, 100, t))/sqrt(nsamples))))
colnames(lag_100_tv_ub) <- c('t', 'means', 'sd')
lag_150_tv_ub <- data.frame(cbind(c(1:200),sapply(1:200, function(t) mean(upper_bound_estimates(meetingtimes_150, 150, t))),
                                  sapply(1:200, function(i) sd(upper_bound_estimates(meetingtimes_150, 150, t))/sqrt(nsamples))))
colnames(lag_150_tv_ub) <- c('t', 'means', 'sd')
lag_500_tv_ub <- data.frame(cbind(c(1:200),sapply(1:200, function(t) mean(upper_bound_estimates(meetingtimes_500, 500, t))),
                                  sapply(1:200, function(i) sd(upper_bound_estimates(meetingtimes_500, 500, t))/sqrt(nsamples))))
colnames(lag_500_tv_ub) <- c('t', 'means', 'sd')


g_ub_PG <- ggplot(data=lag_1_tv_ub,aes(x = t,y=means, col = '1.lag_1'))+geom_line(aes(y=means)) + 
  geom_line(data=lag_10_tv_ub, aes(y=means, col = '2.lag_10')) + 
  geom_line(data=lag_50_tv_ub, aes(y=means, col = '3.lag_50')) +
  geom_line(data=lag_100_tv_ub, aes(y=means, col = '4.lag_100')) + 
  geom_line(data=lag_150_tv_ub, aes(y=means, col = '5.lag_150')) +
  geom_line(data=lag_500_tv_ub, aes(y=means, col = '6.lag_500')) + 
  ylim(0,1) + xlim(0,200) + theme(legend.position="bottom", legend.box = "horizontal") + labs(col="Lag") +
  labs(y = TeX("d_{TV}( pi_t , pi ) \ bounds")) + labs(title='Polya-Gamma')


#mcmc_traceplot <- mcmc_traceplot + xlim(0,200) 
#tv_bound_pic_1 <- grid.arrange(mcmc_traceplot, g_ub, nrow=2)



