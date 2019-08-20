# library(devtools)
# install_github("pierrejacob/debiasedmcmc")

# load packages
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
registerDoParallel(cores = detectCores())



rm(list = ls())
set.seed(21)

source("tv_wasserstein_bounds.R")

# Functions
samplemeetingtime <- function(single_kernel, coupled_kernel, rinit, max_iterations = Inf, L=1){
  
  if(max_iterations != Inf)
  {
    chain1 <- rep(0, max(max_iterations,L))
    chain2 <- rep(0, max(max_iterations,L))
  } else {
    chain1 <- rep(0, 1000)
    chain2 <- rep(0, 1000)
  }
  
  
  chain_state1 <- rinit()
  chain_state2 <- rinit()
  
  chain1[1] <- chain_state1
  chain2[1] <- chain_state2
  
  if (L>1){
    for (l in 1:L)
    {
      chain1[1] <- single_kernel(chain1[1])
    }
  }
  
  chain1[2] <- single_kernel(chain1[1])
  chain_state1 <- chain1[2]
  
  iter <- L
  meet <- FALSE
  finished <- FALSE
  meetingtime <- Inf
  # iter here is L; at this point we have X[2]=X_L,Y[1]=Y_0 and we are going to generate successively X_t,Y_{t-1} where iter = t
  while (!finished && iter < max_iterations){
    iter <- iter + 1
    if (meet){
      chain_state1 <- single_kernel(chain_state1)
      chain_state2 <- chain_state1
    } else {
      res_coupled_kernel <- coupled_kernel(chain_state1, chain_state2)
      chain_state1 <- res_coupled_kernel$chain_state1
      chain_state2 <- res_coupled_kernel$chain_state2
      if (all(chain_state1 == chain_state2) && !meet){
        # recording meeting time tau
        meet <- TRUE
        meetingtime <- iter
      }
    }
    
    chain1[(iter)+2-L] <- chain_state1
    chain2[(iter)+1-L] <- chain_state2
    
    # Appending extra memory for chains
    if ((iter)+2-L+1000>length(chain1))
    {
      chain1 <- c(chain1, rep(0,1000))
      chain2 <- c(chain2, rep(0,1000))
    }
    
    # stop after max(m, tau) steps
    if (iter >= meetingtime){
      finished <- TRUE
    }
  }
  return(list(meetingtime = meetingtime, chain1 = chain1, chain2 = chain2, iteration = iter, finished = finished))
}

target <- function(x){
  evals <- log(0.5) + dnorm(x, mean = c(-4, 4), sd = 1, log = TRUE)
  return(max(evals) + log(sum(exp(evals - max(evals)))))
}

get_pb <- function(sd_proposal, initmean, initsd){
  # Markov kernel of the chain
  single_kernel <- function(chain_state){
    proposal_value <- rnorm(1, mean=chain_state, sd=sd_proposal)
    proposal_pdf <- target(proposal_value)
    current_pdf <- target(chain_state)
    accept <- (log(runif(1)) < (proposal_pdf - current_pdf))
    if (accept){
      return(proposal_value)
    } else {
      return(chain_state)
    }
  }
  
  # Markov kernel of the coupled chain
  coupled_kernel <- function(chain_state1, chain_state2){
    proposal_value <- rnorm_max_coupling(chain_state1, chain_state2, sd_proposal, sd_proposal)
    proposal1 <- proposal_value[1]
    proposal2 <- proposal_value[2]
    proposal_pdf1 <- target(proposal1)
    proposal_pdf2 <- target(proposal2)
    current_pdf1 <- target(chain_state1)
    current_pdf2 <- target(chain_state2)
    logu <- log(runif(1))
    accept1 <- FALSE
    accept2 <- FALSE
    if (is.finite(proposal_pdf1)){
      accept1 <- (logu < (proposal_pdf1 - current_pdf1))
    }
    if (is.finite(proposal_pdf2)){
      accept2 <- (logu < (proposal_pdf2 - current_pdf2))
    }
    if (accept1){
      chain_state1 <- proposal1
    }
    if (accept2){
      chain_state2 <- proposal2
    }
    return(list(chain_state1 = chain_state1, chain_state2 = chain_state2))
  }
  rinit <- function() rnorm(1, initmean, initsd)
  return(list(rinit = rinit, single_kernel = single_kernel, coupled_kernel = coupled_kernel))
}


################################# 
## target distribution
target <- function(x){
  evals <- log(0.5) + dnorm(x, mean = c(-4, 4), sd = 1, log = TRUE)
  return(max(evals) + log(sum(exp(evals - max(evals)))))
}
curve(sapply(x, function(v) exp(target(v))), from = -10, to = 10)

bimodal_pdf_df <- data.frame(sapply(seq(-10,10,length.out = 10000), function(v) exp(target(v))))
colnames(bimodal_pdf_df) <- c('pdf')

################################# Generating data
# Hard setting: bad proposal, bad init
pb <- get_pb(1, initmean = 10, initsd = 1)

# Plotting chain trails
chain_length <- 500
number_of_chains <- 1000

all_chains <- matrix(0, nrow = number_of_chains, ncol = chain_length)
for (i in 1:number_of_chains)
{
  X <- rep(0, chain_length)
  X[1] <- pb$rinit()
  for (t in 2:chain_length){
    X[t] <- pb$single_kernel(X[(t-1)])
  }
  all_chains[i,] <- X
}
# Traceplot data
#write.csv(t(all_chains), "bad_bimodal_traceplots.csv")




## Generating coupled chains
nsamples <- 20000
max_iterations <- Inf

t_start <- 0
t_end <- 10000
meetings_lag_1 <- rep(0, nsamples)
meetings_lag_18000 <- rep(0, nsamples)

for(irep in 1:nsamples){
  lag_1 <- samplemeetingtime(pb$single_kernel, pb$coupled_kernel, pb$rinit, max_iterations,  L = 1)
  meetings_lag_1[irep] <- lag_1$meetingtime
  
  lag_18000 <- samplemeetingtime(pb$single_kernel, pb$coupled_kernel, pb$rinit, max_iterations,  L = 18000)
  meetings_lag_18000[irep] <- lag_18000$meetingtime
  
  print(irep)
}

#write.csv(meetings_lag_1, "bad_bimodal_lag_1.csv")
#write.csv(meetings_lag_18000, "bad_bimodal_lag_18000.csv")



