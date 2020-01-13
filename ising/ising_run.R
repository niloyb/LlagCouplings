library(doParallel)
library(doRNG)

source("coupling.R")
source("ising/ising_functions.R")

set.seed(1)
registerDoParallel(cores = detectCores())

size <- 32
# We're ultimately interested in a low temperature, for which the SSG sampler is not mixing well.
beta <- 0.46
NREP <- 5e2
max_iterations <- 1e7

# Single-site Gibbs.
ssg_lag <- 1e6
ssg_kernel <- get_ising_ssg_kernel(size, beta)
ssg_meeting_times <- foreach(i = 1:NREP, .combine=rbind) %dorng% {
  res <- simulate_meeting_time(ssg_kernel$kernel, ssg_kernel$coupled_kernel, ssg_kernel$init, max_iterations = max_iterations, L = ssg_lag)
  res$meeting_time
}
write.csv(ssg_meeting_times, "ising_ssg_meeting_times.csv")

# Parallel tempering.
pt_lag <- 2e4
# Introduce a grid of lower values of beta.
nchains <- 12
betas_grid <- seq(from = 0.3, to = beta, length.out = nchains)
# Probability of performing a swap move rather than a sweep of single site Gibbs updates.
proba_swapmove <- 1 / 50
pt_cmp <- function(chain_states1, chain_states2) {
  all(sapply(1:length(chain_states1), function(i) all(chain_states1[[i]] == chain_states2[[i]])))
}
pt_kernel <- get_ising_pt_kernel(size, betas_grid, proba_swapmove)
pt_meeting_times <- foreach(i = 1:NREP, .combine=rbind) %dorng% {
  res <- simulate_meeting_time(pt_kernel$kernel, pt_kernel$coupled_kernel, pt_kernel$init, max_iterations = max_iterations, L = pt_lag, cmp = pt_cmp)
  res$meeting_time
}
write.csv(pt_meeting_times, "ising_pt_meeting_times.csv")
