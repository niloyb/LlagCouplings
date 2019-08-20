# Generate TV and Wasserstein bounds given coupling times and the pair of chains

tv_upper_bound_estimates <- function(coupling_times, L, t)
{
  return(pmax(0,ceiling((coupling_times-L-t)/L)))
}

wasserstein_distance_upper_bound_estimate <- function(coupling_time, chain1, chain2, L, t)
{
  J <- ceiling((coupling_time-L-t)/L)
  
  if (J>0){
    #chain1 <- chain1[c(1:(t+2+(J-1)*L))]
    #chain2 <- chain2[c(1:(t+2+(J-1)*L))]
    thinned_times <- t+2 + c(0:(J-1))*L
    w <- sum(abs(chain1[thinned_times]-chain2[(thinned_times-1)]))
  } else {
    w <- 0
  }
  return(w)
}

# Multiple estimates of Wasserstein bounds
wasserstein_distance_upper_bound_estimates <- function(coupling_times, chain1s, chain2s, nsamples, L, t)
{
  repeats <- nsamples
  wasserstein_distance_upper_bound_estimates <- rep(0,repeats)
  for (i in 1: repeats){
    wasserstein_distance_upper_bound_estimates[i] <- wasserstein_distance_upper_bound_estimate(coupling_times[i], chain1s[i,], chain2s[i,], L, t)
  }
  return(wasserstein_distance_upper_bound_estimates)
}

