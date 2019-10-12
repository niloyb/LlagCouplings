library(unbiasedmcmc)

get_ising_ssg_kernel <- function(size, beta) {
  # precomputed probability for single-site flips, given sum of neighbors
  ss <- c(-4, -2, 0, 2, 4)
  probas <- exp(ss * beta) / (exp(ss * beta) + exp(-ss * beta))

  return(list(
    init = function() {
      res <- ising_rinit(size)
      return(list(
        chain_state = res,
        current_pdf = NA
      ))
    },
    kernel = function(chain_state, current_pdf, iteration) {
      res <- ising_single_kernel(chain_state, probas)
      return(list(
        chain_state = res,
        current_pdf = NA
      ))
    },
    coupled_kernel = function(chain_state1, chain_state2, current_pdf1, current_pdf2, iteration) {
      res <- ising_coupled_kernel(chain_state1, chain_state2, probas)
      return(list(
        chain_state1 = res$chain_state1,
        chain_state2 = res$chain_state2,
        current_pdf1 = NA,
        current_pdf2 = NA
      ))
    }
  ))
}

get_ising_pt_kernel <- function(size, betas, proba_swapmove) {
  # The functions here are adapted from ising.R in the unbasedmcmc package to avoid computing
  # unnecessary state which is not used in our results.
  nchains <- length(betas)
  
  # precomputed probability for single-site flips, given sum of neighbors
  ss <- c(-4,-2,0,2,4)
  probas <-  sapply(betas, function(beta) exp(ss * beta) / (exp(ss * beta) + exp(-ss * beta)))

  ising_ssg_kernel <- list()
  for (i in 1:length(betas)) {
    ising_ssg_kernel[[i]] <- get_ising_ssg_kernel(size, betas[i])
  }
    
  ising_pt_rinit <- function() {
    chain_states <- list()
    for (ichain in 1:nchains) {
      res <- ising_ssg_kernel[[ichain]]$init()
      chain_states[[ichain]] <- res$chain_state
    }
    return(list(
      chain_state = chain_states,
      current_pdf = NA
    ))
  }
  
  ising_pt_single_kernel <- function(chain_states, current_pdf, iteration) {
    u_iteration <- runif(1)
    nchains <- length(chain_states)
    if (u_iteration < proba_swapmove) {
      # swap move
      sumstates <- unlist(lapply(chain_states, ising_sum_))
      nswap_attempts <- 1
      for (ichain in 1:(nchains-1)) {
        tXi <- sumstates[ichain]
        tXip1 <- sumstates[ichain+1]
        swapaccept_logprob <- (betas[ichain] - betas[ichain+1]) * (tXip1 - tXi)
        swapaccept_u <- runif(1)
        if (log(swapaccept_u) < swapaccept_logprob) {
          # do swap
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
      for (ichain in 1:nchains) {
        res <- ising_ssg_kernel[[ichain]]$kernel(chain_states[[ichain]], NA, iteration)
        chain_states[[ichain]] <- res$chain_state
      }
    }
    
    return(list(
      chain_state = chain_states,
      current_pdf = NA
    ))
  }
  
  ising_pt_coupled_kernel <- function(chain_states1, chain_states2, current_pdf1, current_pdf2, iteration) {
    nchains <- length(chain_states1)
    u_iteration <- runif(1)
    if (u_iteration < proba_swapmove) {
      # swap move
      sumstates1 <- unlist(lapply(chain_states1, ising_sum_))
      sumstates2 <- unlist(lapply(chain_states2, ising_sum_))
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
          tmp <- chain_states1[[ichain]]
          chain_states1[[ichain]] <- chain_states1[[ichain+1]]
          chain_states1[[ichain+1]] <- tmp
          tmp <- sumstates1[ichain]
          sumstates1[ichain] <- sumstates1[ichain+1]
          sumstates1[ichain+1] <- tmp
        }
        if (log(swapaccept_u) < swapaccept_logprob2){
          # do swap
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
        res <- ising_ssg_kernel[[ichain]]$coupled_kernel(chain_states1[[ichain]], chain_states2[[ichain]], NA, NA, iteration)
        chain_states1[[ichain]] <- res$chain_state1
        chain_states2[[ichain]] <- res$chain_state2
      }
    }
    return(list(
      chain_state1 = chain_states1,
      chain_state2 = chain_states2,
      current_pdf1 = NA,
      current_pdf2 = NA
    ))
  }
  
  return(list(
    init = ising_pt_rinit,
    kernel = ising_pt_single_kernel,
    coupled_kernel = ising_pt_coupled_kernel
  ))
}
