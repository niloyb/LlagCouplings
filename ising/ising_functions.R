
# initialization of Gibbs chain for Ising model
# by flipping a fair coin independently for each site in the size x size grid 
ising_rinit <- function(size = 32){
  initial_values <- 2 * rbinom(size*size, 1, 0.5) - 1
  state <- matrix(initial_values, nrow = size, ncol = size)
  return(state)
}

# one step of Gibbs sampler, i.e. one full sweep over all components
# probas should be a vector of length 5, containing proba of drawing +1
# given that the sum of neighboring spins is {-4,-2,0,+2,+4}
ising_single_kernel <- function(chain_state, probas){
  chain_state <- ising_gibbs_sweep_(chain_state, probas)
  return(chain_state)
}

# one step of coupled Gibbs sampler, i.e. one full sweep over all components
# probas should be a vector of length 5, containing proba of drawing +1
# given that the sum of neighboring spins is {-4,-2,0,+2,+4}
ising_coupled_kernel <- function(chain_state1, chain_state2, probas){
  res_ <- ising_coupled_gibbs_sweep_(chain_state1, chain_state2, probas)
  return(list(chain_state1 = res_$state1, chain_state2 = res_$state2))
}

# initialization of parallel tempering Gibbs chain for Ising model
ising_pt_rinit <- function(nchains, size = 32){
  chain_states <- list()
  for (ichain in 1:nchains){
    chain_states[[ichain]] <- ising_rinit(size)
  }
  return(chain_states)
}

# one step of parallel tempering Gibbs
# chain_states should be a list containing a number N of chains
# sumstates should be a vector containing N integers corresponding to the sum of spins as computed by ising_sum_
# betas should be a vector of N inverse temperatures
# probs should be a 5 x N matrix of probabilities of drawing +1 given that the sum of neighboring spins is {-4,-2,0,+2,+4}
# proba_swapmove is a number in (0,1) representing the probability of performing swap moves
ising_pt_single_kernel <- function(chain_states, sumstates, betas, probas, proba_swapmove){
  u_iteration <- runif(1)
  nchains <- length(chain_states)
  nswap_attempts <- 0
  nswap_accepts <- rep(0, nchains-1)
  if (u_iteration < proba_swapmove){
    # swap move
    nswap_attempts <- 1
    for (ichain in 1:(nchains-1)){
      tXi <- sumstates[ichain]
      tXip1 <- sumstates[ichain+1]
      swapaccept_logprob <- (betas[ichain] - betas[ichain+1]) * (tXip1 - tXi)
      swapaccept_u <- runif(1)
      if (log(swapaccept_u) < swapaccept_logprob){
        # do swap
        nswap_accepts[ichain] <- 1
        tmp <- chain_states[[ichain]]
        chain_states[[ichain]] <- chain_states[[ichain+1]]
        chain_states[[ichain+1]] <- tmp
        tmp <- sumstates[ichain]
        sumstates[ichain] <- sumstates[ichain+1]
        sumstates[ichain+1] <- tmp
      }
    }
  } else {
    # Gibbs move
    for (ichain in 1:nchains){
      chain_states[[ichain]] <- ising_gibbs_sweep_(chain_states[[ichain]], proba_beta = probas[,ichain])
    }
  }
  sumstates <- unlist(lapply(chain_states, ising_sum_))
  return(list(chain_states = chain_states, sumstates = sumstates, nswap_attempts = nswap_attempts, nswap_accepts = nswap_accepts))
}


# one step of coupled parallel tempering Gibbs
# chain_states1,chain_states2 should be lists containing a number N of chains
# sumstates1,sumstates2 should be vectors containing N integers corresponding to the sum of spins as computed by ising_sum_
# betas should be a vector of N inverse temperatures
# probs should be a 5 x N matrix of probabilities of drawing +1 given that the sum of neighboring spins is {-4,-2,0,+2,+4}
# proba_swapmove is a number in (0,1) representing the probability of performing swap moves
ising_pt_coupled_kernel <- function(chain_states1, chain_states2, sumstates1, sumstates2, betas, probas, proba_swapmove){
  nchains <- length(chain_states1)
  nswap_attempts <- 0
  nswap_accepts1 <- rep(0, nchains-1)
  nswap_accepts2 <- rep(0, nchains-1)
  u_iteration <- runif(1)
  if (u_iteration < proba_swapmove){
    # swap move
    nswap_attempts <- 1
    for (ichain in 1:(nchains-1)){
      tXi_1 <- sumstates1[ichain]
      tXip1_1 <- sumstates1[ichain+1]
      tXi_2 <- sumstates2[ichain]
      tXip1_2 <- sumstates2[ichain+1]
      deltabeta <- betas[ichain] - betas[ichain+1]
      # swapaccept_logprob <- tXi * (-deltabeta) + tXip1 * deltabeta
      swapaccept_logprob1 <- deltabeta * (tXip1_1 - tXi_1)
      swapaccept_logprob2 <- deltabeta * (tXip1_2 - tXi_2)
      swapaccept_u <- runif(1)
      if (log(swapaccept_u) < swapaccept_logprob1){
        # do swap
        nswap_accepts1[ichain] <- 1
        tmp <- chain_states1[[ichain]]
        chain_states1[[ichain]] <- chain_states1[[ichain+1]]
        chain_states1[[ichain+1]] <- tmp
        tmp <- sumstates1[ichain]
        sumstates1[ichain] <- sumstates1[ichain+1]
        sumstates1[ichain+1] <- tmp
      }
      if (log(swapaccept_u) < swapaccept_logprob2){
        # do swap
        nswap_accepts2[ichain] <- 1
        tmp <- chain_states2[[ichain]]
        chain_states2[[ichain]] <- chain_states2[[ichain+1]]
        chain_states2[[ichain+1]] <- tmp
        tmp <- sumstates2[ichain]
        sumstates2[ichain] <- sumstates2[ichain+1]
        sumstates2[ichain+1] <- tmp
      }
    }
  } else {
    # Gibbs move
    for (ichain in 1:nchains){
      res_ <- ising_coupled_gibbs_sweep_(chain_states1[[ichain]], chain_states2[[ichain]], probas[,ichain])
      chain_states1[[ichain]] <- res_$state1
      chain_states2[[ichain]] <- res_$state2
    }
  }
  sumstates1 <- unlist(lapply(chain_states1, ising_sum_))
  sumstates2 <- unlist(lapply(chain_states2, ising_sum_))
  return(list(chain_states1 = chain_states1, chain_states2 = chain_states2,
              sumstates1 = sumstates1, sumstates2 = sumstates2,
              nswap_attempts = nswap_attempts,
              nswap_accepts1 = nswap_accepts1, nswap_accepts2 = nswap_accepts2))
}

ising_gibbs_meetingtime <- function(beta, size = 32, lag = 1, max_iterations = Inf){
  ss_ <- c(-4,-2,0,2,4)
  # precomputed probability for single-site flips, given sum of neighbors
  proba_ <-  exp(ss_*beta) / (exp(ss_*beta) + exp(-ss_*beta))
  # initialize
  chain_state1 <- ising_rinit(size)
  chain_state2 <- ising_rinit(size)
  iter <- 0 # (index of current time for chain1)
  # move first chain
  for (t in 1:lag){
    iter <- iter + 1
    chain_state1 <- ising_single_kernel(chain_state1, proba_)
  }
  # iterate
  meetingtime <- Inf
  # iter here is 1; at this point we have X_1,Y_0 and we are going to generate successively X_t,Y_{t-1} where iter = t
  while (is.infinite(meetingtime) && iter < max_iterations){
    iter <- iter + 1
    # use coupled kernel
    res_ <- ising_coupled_kernel(chain_state1, chain_state2, proba_)
    chain_state1 <- res_$chain_state1
    chain_state2 <- res_$chain_state2
    # check if meeting happens
    allequal <- all(chain_state1 == chain_state2)
    if (allequal){
      # recording meeting time tau
      meetingtime <- iter
    }
  }
  # return meeting time
  return(list(meetingtime = meetingtime, iteration = iter))
}


# The function below draws meeting times 
# using "coupled parallel tempering" for the inverse temperatures provided in the vector "betas"
ising_pt_meetingtime <- function(betas, proba_swapmove, size = 32, lag = 1, max_iterations = Inf){
  nchains  <- length(betas)
  ss_ <- c(-4,-2,0,2,4)
  probas_ <-  sapply(betas, function(beta) exp(ss_*beta) / (exp(ss_*beta) + exp(-ss_*beta)))
  # initialize
  chain_states1 <- ising_pt_rinit(nchains, size)
  chain_states2 <- ising_pt_rinit(nchains, size)
  sumstates1 <- unlist(lapply(chain_states1, ising_sum_))
  sumstates2 <- unlist(lapply(chain_states2, ising_sum_))
  # move first chain
  iter <- 1
  for (t in 1:lag){
    iter <- iter + 1
    # chain_state1 <- ising_single_kernel(chain_state1, proba_)
    res_single_kernel <- ising_pt_single_kernel(chain_states1, sumstates1, betas, probas_, proba_swapmove)
    chain_states1 <- res_single_kernel$chain_states
    sumstates1 <- res_single_kernel$sumstates
  }
  # iterate
  meetingtime <- Inf
  while (is.infinite(meetingtime) && iter < max_iterations){
    iter <- iter + 1
    # use coupled kernel
    res_ <- ising_pt_coupled_kernel(chain_states1, chain_states2, sumstates1, sumstates2, betas, probas_, proba_swapmove)
    chain_states1 <- res_$chain_states1
    chain_states2 <- res_$chain_states2
    sumstates1 <- res_$sumstates1
    sumstates2 <- res_$sumstates2
    # check if meeting happens
    allequal <- all(sapply(1:nchains, function(i) all(chain_states1[[i]] == chain_states2[[i]])))
    if (allequal){
      # recording meeting time tau
      meetingtime <- iter
    }
  }
  # return meeting time
  return(list(meetingtime = meetingtime, iteration = iter))
}
