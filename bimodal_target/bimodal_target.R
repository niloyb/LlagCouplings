library(doParallel)
library(doRNG)
#library(unbiasedmcmc)

source("coupling.R")

set.seed(1)
registerDoParallel(cores = detectCores())

# Log density of the bimodal target: 0.5*N(-4,1)+0.5*N(+4,1)
logtarget <- function(x) {
  evals <- log(0.5) + dnorm(x, mean = c(-4, 4), sd = 1, log = TRUE)
  return(max(evals) + log(sum(exp(evals - max(evals)))))
}

# Sample from a maximal coupling of two distributions p and q
# Distributions p and q are specified by random number generators and log-pdfs
get_max_coupling <- function(rp, dp, rq, dq){
  function(){
    x <- rp(1)
    if (dp(x) + log(runif(1)) < dq(x)){
      return(list(xy = c(x,x), identical = TRUE))
    } else {
      reject <- TRUE
      y <- NA
      while (reject){
        y <- rq(1)
        reject <- (dq(y) + log(runif(1)) < dp(y))
      }
      return(list(xy = c(x,y), identical = FALSE))
    }
  }
}

# Samples from maximal coupling of two univariate Normal distributions
# Normal distributions are specified through their means and standard deviations
rnorm_max_coupling <- function(mu1, mu2, sigma1, sigma2){
  f <- get_max_coupling(function(n) rnorm(n, mu1, sigma1),
                        function(x) dnorm(x, mu1, sigma1, log = TRUE),
                        function(n) rnorm(n, mu2, sigma2),
                        function(x) dnorm(x, mu2, sigma2, log = TRUE))
  return(f())
}

# get_pb returns the kernels of a single RWMH and a coupled RWMH chain
get_pb <- function(sd_proposal, initmean, initsd) {
  # Initilization the RWMH chain
  rinit <- function() {
    chain_state <- rnorm(1, initmean, initsd)
    return(list(
      chain_state = chain_state,
      current_pdf = logtarget(chain_state)
    ))
  }
  
  # Markov kernel of the chain
  single_kernel <- function(chain_state, current_pdf, iteration) {
    proposal_state <- rnorm(1, mean=chain_state, sd=sd_proposal)
    proposal_pdf <- logtarget(proposal_state)
    if (log(runif(1)) < (proposal_pdf - current_pdf)) {
      chain_state <- proposal_state
      current_pdf <- proposal_pdf
    }
    return(list(
      chain_state = chain_state,
      current_pdf = current_pdf
    ))
  }
  
  # Markov kernel of the coupled chain
  coupled_kernel <- function(chain_state1, chain_state2, current_pdf1, current_pdf2, iteration) {
    proposal_res <- rnorm_max_coupling(chain_state1, chain_state2, sd_proposal, sd_proposal)
    proposal1 <- proposal_res$xy[1]
    proposal2 <- proposal_res$xy[2]
    proposal_pdf1 <- logtarget(proposal1)
    proposal_pdf2 <- logtarget(proposal2)
    logu <- log(runif(1))
    if (is.finite(proposal_pdf1)){
      if (logu < (proposal_pdf1 - current_pdf1)) {
        chain_state1 <- proposal1
        current_pdf1 <- proposal_pdf1
      }
    }
    if (is.finite(proposal_pdf2)){
      if (logu < (proposal_pdf2 - current_pdf2)) {
        chain_state2 <- proposal2
        current_pdf2 <- proposal_pdf2
      }
    }
    return(list(
      chain_state1 = chain_state1, 
      chain_state2 = chain_state2,
      current_pdf1 = current_pdf1,
      current_pdf2 = current_pdf2
    ))
  }
 
  return(list(
    rinit = rinit,
    single_kernel = single_kernel,
    coupled_kernel = coupled_kernel
  ))
}


# Create kernels for chain beginning far from target distribution.
pb <- get_pb(1, initmean = 10, initsd = 1)

# Generate MCMC chains.
chain_length <- 500
number_of_chains <- 1000

# Single RWMH chain
all_chains <- foreach(i = 1:number_of_chains, .combine = rbind) %dorng% {
    res <- pb$rinit()
    x <- res$chain_state
    for (t in 2:chain_length) {
      res <- pb$single_kernel(res$chain_state, res$current_pdf, t)
      x <- res$chain_state
    }
    x
  }
write.csv(all_chains, "bad_bimodal_500.csv")

# Generating coupled chain meeting times.
nsamples <- 10000
max_iterations <- Inf

# Meetings times are under the choice of a small (1) and large (18000) lag
for (lag in c(1, 18000)) {
  meeting_times <- foreach(i = 1:nsamples, .combine = rbind) %dorng% {
    res <- simulate_meeting_time(pb$single_kernel, pb$coupled_kernel, pb$rinit, max_iterations = max_iterations, L = lag)
    res$meeting_time
  }
  write.csv(meeting_times, sprintf("bad_bimodal_lag=%d.csv", lag))
}
