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
library(doRNG)
library(doParallel)

source("tv_wasserstein_bounds.R")


rm(list = ls())
set.seed(1)
registerDoParallel(cores = detectCores())




#######################################################################   Markov Chain TV bounds plot: simple MH targetting N(0,1)
chain_length <- 210
number_of_chains <- 50000
mh_std <- 0.5
X <- rep(0, chain_length)
X[1] <- 10

all_chains <- matrix(0, nrow = number_of_chains, ncol = chain_length)
for (i in 1:number_of_chains)
{
  for (t in 2:chain_length){
    proposal <- rnorm(1, X[(t-1)], sd = mh_std)
    X[t] <- X[(t-1)]
    if (log(runif(1)) < - 0.5*(proposal^2-X[(t-1)]^2)) { X[t] <- proposal }
  }
  
  all_chains[i,] <- X
}

# Traceplot
# matplot(t(all_chains), type = 'l')


# ggridges Traceplot
trace_plot_df <- data.frame(cbind(c(1:chain_length),t(all_chains)))
trace_plot_df <- melt(trace_plot_df, id="X1")
colnames(trace_plot_df) <- c('t', 'chain', 'value')
#g_1 <- ggplot(data=trace_plot_df,aes(x=t, y=X, colour=chain)) + geom_line() + theme(legend.position = "none")
#mcmc_traceplot <- g_1
g_1 <- ggplot(data=trace_plot_df,aes(x=value, y = factor(t))) + geom_density_ridges(scale = 10)
g_1 <- g_1 + ylim(-50, 250)
g_1 <- g_1 + scale_x_continuous() + scale_y_discrete(breaks = c(-50, 0,50,100,150, 200, 250))
g_1 <- g_1 + xlab("x") + ylab("t") + xlim(-5, 11)
g_1 <- g_1 + coord_flip()
mcmc_traceplot <- g_1
mcmc_traceplot








#######################################################################   Coupling algorithm
logtarget <- function(x) dnorm(x, mean = 0, sd = 1, log = TRUE)
rnorm_max_coupling <- function(mu1, mu2, sigma1, sigma2){
  x <- rnorm(1, mu1, sigma1)
  if (dnorm(x, mu1, sigma1, log = TRUE) + log(runif(1)) < dnorm(x, mu2, sigma2, log = TRUE)){
    return(c(x,x))
  } else {
    reject <- TRUE
    y <- NA
    while (reject){
      y <- rnorm(1, mu2, sigma2)
      reject <- (dnorm(y, mu2, sigma2, log = TRUE) + log(runif(1)) < dnorm(y, mu1, sigma1, log = TRUE))
    }
    return(c(x,y))
  }
}
coupledMH_kernel <- function(chain_state1, chain_state2, sd_proposal){
  proposal_value <- rnorm_max_coupling(chain_state1, chain_state2, sd_proposal, sd_proposal)
  proposal_pdf1 <- logtarget(proposal_value[1])
  proposal_pdf2 <- logtarget(proposal_value[2])
  current_pdf1 <- logtarget(chain_state1)
  current_pdf2 <- logtarget(chain_state2)
  logu <- log(runif(1))
  if (is.finite(proposal_pdf1)){
    if (logu < (proposal_pdf1 - current_pdf1)){
      chain_state1 <- proposal_value[1]
    }
  }
  if (is.finite(proposal_pdf2)){
    if(logu < (proposal_pdf2 - current_pdf2)){
      chain_state2 <- proposal_value[2]
    }
  }
  return(list(chain_state1 = chain_state1, chain_state2 = chain_state2))
}
meeting_time <- function(initial_value, mh_std, lag = 1, max_iterations = 1000)
{
  chain1 <- rep(0, max_iterations)
  chain2 <- rep(0, max_iterations)
  
  chain1[1] <- initial_value
  chain2[1] <- initial_value
  
  # Initialising such that X[1] is the (L-1)^th marginal
  if(lag>1)
  {
    for (l in 1:(lag-1))
    {
      chain1_proposal <- rnorm(1, chain1[1], sd = mh_std)
      if (log(runif(1)) < - 0.5*(chain1_proposal^2-chain1[1]^2)) { chain1[1] <- chain1_proposal }
    }
  }
  
  
  chain1[2] <- chain1[1]
  chain1_proposal <- rnorm(1, chain1[1], sd = mh_std)
  if (log(runif(1)) < - 0.5*(chain1_proposal^2-chain1[1]^2)) { chain1[2] <- chain1_proposal }
  
  not_coupled <- TRUE
  tau <- NA
  t <- 2
  while(not_coupled && t<max_iterations)
  {
    sample <- coupledMH_kernel(chain1[t], chain2[(t-1)], mh_std)
    chain1[(t+1)] <- sample$chain_state1
    chain2[t] <- sample$chain_state2
    if (chain1[(t+1)] == chain2[t])
    {
      tau <- t - 1 + lag
      not_coupled <- FALSE
    }
    t <- t+1
  } 
  
  return(list('meeting_time'=tau, 'chain1'=chain1, 'chain2'=chain2))
}


###### We adaptively choose L values based on Algo in paper
repeats <- 10000
max_iterations <- 1000

meeting_1_lag <- rep(0, repeats)
chain1_1_lag <- matrix(0, nrow = repeats, ncol = max_iterations)
chain2_1_lag <- matrix(0, nrow = repeats, ncol = max_iterations)
meeting_150_lag <- rep(0, repeats)
chain1_150_lag <- matrix(0, nrow = repeats, ncol = max_iterations)
chain2_150_lag <- matrix(0, nrow = repeats, ncol = max_iterations)

for (i in 1:repeats){
  lag_1 <- meeting_time(10, mh_std=0.5, lag = 1, max_iterations)
  meeting_1_lag[i] <- lag_1$meeting_time
  chain1_1_lag[i,] <- lag_1$chain1
  chain2_1_lag[i,] <- lag_1$chain2
  
  lag_150 <- meeting_time(10, mh_std=0.5, lag = 150, max_iterations)
  meeting_150_lag[i] <- lag_150$meeting_time
  chain1_150_lag[i,] <- lag_150$chain1
  chain2_150_lag[i,] <- lag_150$chain2
  
  print(i)
}

######################## ggplot for TV and Wasserstein bounds
lag_1_ub <- data.frame(cbind((c(1:chain_length)-1), sapply((c(1:chain_length)-1), function(t) mean(tv_upper_bound_estimates(meeting_1_lag, 1, t))),
                           sapply((c(1:chain_length)-1), function(t) mean(wasserstein_distance_upper_bound_estimates(meeting_1_lag, chain1_1_lag, chain2_1_lag, repeats, 1, t)))))
colnames(lag_1_ub) <- c('t', 'TV', 'Wass')

lag_150_ub <- data.frame(cbind((c(1:chain_length)-1), sapply((c(1:chain_length)-1), function(t) mean(tv_upper_bound_estimates(meeting_150_lag, 150, t))),
                             sapply((c(1:chain_length)-1), function(t) mean(wasserstein_distance_upper_bound_estimates(meeting_150_lag, chain1_150_lag, chain2_150_lag,repeats, 150, t)))))
colnames(lag_150_ub) <- c('t', 'TV', 'Wass')

######################## Exact TV calculation from kernel density estimates
nmcmc <- chain_length
X <- t(all_chains)
target <- function(x){
  evals <- dnorm(x, mean = 0, sd = 1, log = TRUE)
  return(max(evals) + log(sum(exp(evals - max(evals)))))
}

## compute TV at each step, using kernel density estimates
## if some numerical error arises, then the TV is set to zero
tvs <- rep(0, chain_length)
for(imcmc in 1:chain_length){
  # compute kernel density estimate, based on Markov chains
  kde_ <- density(X[imcmc,])
  # compute |f(x) - pi(x)| where f is the kernel density estimate based on samples, and pi is the target
  diff_12 <- function(x){
    sapply(x, function(v){
      if ((v > min(kde_$x)) && (v < max(kde_$x))){
        return(abs(approx(kde_$x, kde_$y, xout = v)$y - exp(target(v))))
      } else {
        return(abs(exp(target(v))))
      }
    })
  }
  # plot curve of |f(x) - pi(x)| against x
  curve(diff_12(x), from = -5, to = 15, log = "y", main = paste("iteration", imcmc))
  # try to perform numerical integration
  result_ <- try(integrate(diff_12, lower = -5, upper = 15, stop.on.error = FALSE)$value / 2)
  if (inherits(result_, "try-error")){
    result_ <- 0
  }
  tvs[imcmc] <- result_
}

# plot TV against iterations
plot(x = (c(1:chain_length)-1), y = tvs, type = "l")
# plot TV against iterations, log-scale for the y axis
plot(x = (c(1:chain_length)-1), y = tvs, type = "l", log = "y")

######################## Exact Wasserstein calculations
exact_wasserstein <- rep(0, chain_length)

for (t in 1:chain_length){
  std_norm_samples <- rnorm(number_of_chains)
  std_norm_samples <- sort(std_norm_samples)
  marginal_samples <- all_chains[,t]
  marginal_samples <- sort(marginal_samples)
  
  exact_wasserstein[t] <- mean(abs(std_norm_samples-marginal_samples))
    
}

# plot Wasserstein against iterations
plot(x = (c(1:chain_length)-1), y = exact_wasserstein, type = "l")
# plot Wasserstein against iterations, log-scale for the y axis
plot(x = (c(1:chain_length)-1), y = exact_wasserstein, type = "l", log = "y")


######################## Final plots with MCMC traces, TV bounds and Exact TV calculation from kernel density estimates
exact_TV <- data.frame(cbind((c(1:chain_length)-1), tvs))
g_tv_ub <- ggplot(data=lag_1_ub,aes(x = t,y=TV, linetype = '1.L=1'))+geom_line(aes(y=TV)) + 
  geom_line(data=lag_150_ub, aes(y=TV, linetype = '2.L=150')) + 
  geom_line(data=exact_TV, aes(x= V1, y=tvs, linetype = '3.Exact TV')) + 
  scale_x_continuous(expand = c(0, 0), limits = c(0,210)) + ylim(0,1.5) +
  labs(y = TeX("d_{TV}")) + labs(x = "iteration") + labs(linetype="Lag") +
  scale_linetype_manual(values = c(3:1), breaks=c("1.L=1", "2.L=150", "3.Exact TV"), labels=c("L=1", "L=150", "Exact TV")) +
  theme_grey(base_size = 14) + 
  theme(legend.position="bottom", legend.box = "horizontal") 
  

exact_Wass <- data.frame(cbind((c(1:chain_length)-1), exact_wasserstein))
g_wass_ub <- ggplot(data=lag_1_ub,aes(x = t,y=Wass, linetype = '1.L=1'))+geom_line(aes(y=Wass)) + 
  geom_line(data=lag_150_ub, aes(y=Wass, linetype = '2.L=150')) + 
  geom_line(data=exact_Wass, aes(x= V1, y=exact_wasserstein, linetype = '3.Exact Wasserstein')) + 
  scale_x_continuous(expand = c(0, 0), limits = c(0,210)) + ylim(0,12) +
  labs(y = TeX("d_{W}")) + labs(x = "iteration") + labs(linetype="Lag") +
  scale_linetype_manual(values = c(3:1), breaks=c("1.L=1", "2.L=150", "3.Exact Wasserstein"), labels=c("L=1", "L=150", "Exact Wasserstein")) +
  theme_grey(base_size = 14) + 
  theme(legend.position="bottom", legend.box = "horizontal") 
  


tv_bound_pic_1 <- grid.arrange(mcmc_traceplot, g_tv_ub, g_wass_ub, nrow=3)

#ggsave(filename = "bound_standard_normal.pdf", plot = tv_bound_pic_1, width = 8, height = 6)



