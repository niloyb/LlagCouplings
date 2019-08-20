### Ising model 
rm(list = ls())
set.seed(21)
#
library(doRNG)
library(doParallel)
registerDoParallel(cores = detectCores()-1)
library(Rcpp)
sourceCpp("ising.cpp")
source("ising_functions.R")

## now obtain TV bounds for this Gibbs sampler
source("tv_wasserstein_bounds.R")
## meeting time of single site Gibbs sampler (SSG)
get_ssg_meetingtimes <- function(beta, size, lag, NREP, max_iterations = 1e7){
  ssg_meetings <- foreach(irep = 1:NREP) %dorng% {
    ising_gibbs_meetingtime(beta, size, lag = lag, max_iterations = max_iterations)
  }
  meetingtimes <- sapply(ssg_meetings, function(x) x$meetingtime)
  return(meetingtimes)
}


## we can compare the performance of single site Gibbs to parallel tempering with single site Gibbs moves
# size of the grid
size <- 32
# suppose we're ultimately interested in a low temperature, for which the SSG sampler is not mixing well
beta <- 0.46
# introduce a grid of lower values of beta:
nchains <- 12
betas_grid <- seq(from = 0.3, to = beta, length.out = nchains)
# probability of performing a swap move rather than a sweep of single site Gibbs updates 
proba_swapmove <- 1/50

## we can obtain reasonable meeting times using the parallel tempering strategy
get_pt_meetingtimes <- function(betas_grid, proba_swapmove, size, lag, NREP, max_iterations = 1e7){
  pt_meetings <- foreach(irep = 1:NREP) %dorng% {
    ising_pt_meetingtime(betas_grid, proba_swapmove, size, lag = lag, max_iterations = max_iterations)
  }
  pt_meetingtimes <- sapply(pt_meetings, function(x) x$meetingtime)
  return(pt_meetingtimes)
}

# ## Uncomment the following to generate the file "ising.tvbounds.pt.vs.ssg.RData"
NREP <- 5e2
# ## pt_meetingtimes <- get_pt_meetingtimes(betas_grid, proba_swapmove, size, lag = 1, NREP)
# ## summary(pt_meetingtimes)
pt_lag <- 2e4
pt_meetingtimes2 <- get_pt_meetingtimes(betas_grid, proba_swapmove, size, lag = pt_lag, NREP)
summary(pt_meetingtimes2)
plot(sapply(1:5e4, function(x) mean(tv_upper_bound_estimates(pt_meetingtimes2, pt_lag, t = x))), type = "l", xlab = "iteration", ylab = "tv bounds")
save(NREP, beta, betas_grid, proba_swapmove, nchains, pt_meetingtimes2, pt_lag, file = "ising.tvbounds.pt.RData")
## now to compare, get TV bounds for SSG
## ssg_meetingtimes <- get_ssg_meetingtimes(beta, size, lag = 1, NREP = NREP, max_iterations = 5e7)
## summary(ssg_meetingtimes)
ssg_lag <- 1e6
ssg_meetingtimes2 <- get_ssg_meetingtimes(beta, size, lag = ssg_lag, NREP = NREP, max_iterations = 1e8)

plot(sapply(1:1e6, function(x) mean(tv_upper_bound_estimates(ssg_meetingtimes2, ssg_lag, t = x))), type = "l", xlab = "iteration", ylab = "tv bounds", log = "x")
save(NREP, beta, ssg_meetingtimes2, ssg_lag, file = "ising.tvbounds.ssg.RData")
