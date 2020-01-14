library(doParallel)
library(doRNG)
#library(unbiasedmcmc)


source("coupling.R")
source("ula_versus_mala/helper_functions/mvnorm.R")
source("ula_versus_mala/helper_functions/mvnorm_couplings.R")
Rcpp::sourceCpp("ula_versus_mala/helper_functions/mvnorm.cpp")

set.seed(1)
registerDoParallel(cores = detectCores())

run_for_dimension <- function(dimension) {
  # define covariance matrix of target
  Sigma_pi <- diag(1, dimension, dimension)
  alpha <- 0.5
  for (i in 1:dimension) {
    for (j in 1:dimension) {
      Sigma_pi[i, j] <- alpha ^ (abs(i - j))
    }
  }
  
  Sigma_pi_chol <- chol(Sigma_pi)
  inv_Sigma_pi_chol <- solve(Sigma_pi_chol)
  inv_Sigma_pi <- solve(Sigma_pi)
  
  # target is Multivariate Normal
  # target log-density
  target <- function(x) {
    fast_dmvnorm_chol_inverse(matrix(x, nrow = 1), rep(0, dimension), inv_Sigma_pi_chol)
  }
  
  # gradient of target log-density
  gradtarget <- function(x) {
    -(inv_Sigma_pi %*% x)[,1]
  }
  
  # proposal covariance matrix
  Sigma_proposal <- diag(1, dimension, dimension)
  Sigma_chol <- chol(Sigma_proposal)
  inv_Sigma_chol <- solve(Sigma_chol)
  
  
  # Initial distribution.
  rinit <- function() {
    state <- fast_rmvnorm(1, rep(1, dimension), diag(dimension))[1,]
    return(list(
      chain_state = state,
      current_pdf = target(state)
    ))
  }
  
  #####################################
  # Unadjusted Langevin Algorithm (ULA)
  #####################################
  ula_stepsize <- 0.1 * (1 / dimension)^(1/6)
  ula_single_kernel <- function(chain_state, current_pdf, iteraion) {
    mean_proposal <- chain_state + ula_stepsize / 2 * t(Sigma_chol %*% gradtarget(chain_state))[1,]
    proposal_value <- fast_rmvnorm_chol(1, mean_proposal, sqrt(ula_stepsize) * Sigma_chol)[1,]
    return(list(
      chain_state = proposal_value,
      current_pdf = NA
    ))
  }
  
  ula_coupled_kernel <- function(chain_state1, chain_state2, current_pdf1, current_pdf2, iteration) {
    mean_proposal1 <- chain_state1 + ula_stepsize/2 * t(Sigma_chol %*% gradtarget(chain_state1))[1,]
    mean_proposal2 <- chain_state2 + ula_stepsize/2 * t(Sigma_chol %*% gradtarget(chain_state2))[1,]
    proposal_value <- rmvnorm_reflectionmax(mean_proposal1, mean_proposal2, sqrt(ula_stepsize) * Sigma_chol, inv_Sigma_chol / sqrt(ula_stepsize))
    proposal1 <- proposal_value$xy[,1]
    proposal2 <- proposal_value$xy[,2]
    return(list(
      chain_state1 = proposal1,
      chain_state2 = proposal2,
      current_pdf1 = NA,
      current_pdf2 = NA
    ))
  }

  # Time a single chain.
  niterations <- 1000
  res <- rinit()
  chain_state <- res$chain_state
  current_pdf <- target(chain_state)
  start_time <- proc.time()
  for (iteration in 1:niterations) {
    res <- ula_single_kernel(chain_state, current_pdf, iteration)
    chain_state <- res$chain_state
    current_pdf <- res$current_pdf
  }
  ula_iteration_time <- as.numeric((proc.time() - start_time)[3]) / niterations
  write.csv(ula_iteration_time, sprintf("ula_time_d=%d.csv", dimension))
  
  # Simulate meeting times.
  ula_lag <- 10000
  nsamples <- 50
  max_iterations <- 100000
  ula_meetings <- foreach(irep = 1:nsamples, .combine = rbind) %dorng% {
    res <- simulate_meeting_time(ula_single_kernel, ula_coupled_kernel, rinit, max_iterations = max_iterations, L = ula_lag)
    res$meeting_time
  }
  write.csv(ula_meetings, sprintf("ula_meetings_d=%d_lag=%d.csv", dimension, ula_lag))
  
  ###############################################
  # Metropolis adjusted Langevin Algorithm (MALA)
  ###############################################
  mala_stepsize = (1/dimension)^(1/6)
  mala_single_kernel <- function(chain_state, current_pdf, stepsize) {
    mean_proposal <- chain_state + mala_stepsize/2 * t(Sigma_chol %*%  gradtarget(chain_state))[1,]
    proposal_value <- fast_rmvnorm_chol(1, mean_proposal, sqrt(mala_stepsize) * Sigma_chol)[1,]
    mean_revert <- proposal_value + mala_stepsize/2  * t(Sigma_chol %*% gradtarget(proposal_value))[1,]
    proposal_pdf <- target(proposal_value)
    accept <- (log(runif(1)) < (proposal_pdf - current_pdf + fast_dmvnorm_chol_inverse(matrix(chain_state, nrow = 1), mean_revert, inv_Sigma_chol / sqrt(mala_stepsize)) -
                                  fast_dmvnorm_chol_inverse(matrix(proposal_value, nrow = 1), mean_proposal, inv_Sigma_chol / sqrt(mala_stepsize))))
    if (accept){
      return(list(chain_state = proposal_value, current_pdf = proposal_pdf))
    } else {
      return(list(chain_state = chain_state, current_pdf = current_pdf))
    }
  }
  
  mala_coupled_kernel <- function(chain_state1, chain_state2, current_pdf1, current_pdf2, stepsize) {
    mean_proposal1 <- chain_state1 + mala_stepsize/2 * t(Sigma_chol %*% gradtarget(chain_state1))[1,]
    mean_proposal2 <- chain_state2 + mala_stepsize/2 * t(Sigma_chol %*% gradtarget(chain_state2))[1,]
    proposal_value <- rmvnorm_reflectionmax(mean_proposal1, mean_proposal2, sqrt(mala_stepsize) * Sigma_chol, inv_Sigma_chol / sqrt(mala_stepsize))
    proposal1 <- proposal_value$xy[,1]
    proposal2 <- proposal_value$xy[,2]
    proposal_pdf1 <- target(proposal1)
    proposal_pdf2 <- target(proposal2)
    
    mean_revert1 <- proposal1 + mala_stepsize/2  * t(Sigma_chol %*% gradtarget(proposal1))[1,]
    mean_revert2 <- proposal2 + mala_stepsize/2  * t(Sigma_chol %*% gradtarget(proposal2))[1,]
    
    logu <- log(runif(1))
    
    accept1 <- FALSE
    accept2 <- FALSE
    if (is.finite(proposal_pdf1)){
      accept1 <- (logu < (proposal_pdf1 - current_pdf1 + fast_dmvnorm_chol_inverse(matrix(chain_state1, nrow = 1), mean_revert1, inv_Sigma_chol / sqrt(mala_stepsize)) -
                            fast_dmvnorm_chol_inverse(matrix(proposal1, nrow = 1), mean_proposal1, inv_Sigma_chol / sqrt(mala_stepsize))))
    }
    if (is.finite(proposal_pdf2)){
      accept2 <- (logu < (proposal_pdf2 - current_pdf2 + fast_dmvnorm_chol_inverse(matrix(chain_state2, nrow = 1), mean_revert2, inv_Sigma_chol / sqrt(mala_stepsize)) -
                            fast_dmvnorm_chol_inverse(matrix(proposal2, nrow = 1), mean_proposal2, inv_Sigma_chol / sqrt(mala_stepsize))))
    }
    if (accept1){
      chain_state1 <- proposal1
      current_pdf1 <- proposal_pdf1
    }
    if (accept2){
      chain_state2 <- proposal2
      current_pdf2 <- proposal_pdf2
    }
    return(list(
      chain_state1 = chain_state1, 
      chain_state2 = chain_state2, 
      current_pdf1 = current_pdf1, 
      current_pdf2 = current_pdf2
    ))
  }
  
  niterations <- 1000
  res <- rinit()
  chain_state <- res$chain_state
  current_pdf <- target(chain_state)
  start_time <- proc.time()
  for (iteration in 1:niterations) {
    res <- mala_single_kernel(chain_state, current_pdf, iteration)
    chain_state <- res$chain_state
    current_pdf <- res$current_pdf
  }
  mala_iteration_time <- as.numeric((proc.time() - start_time)[3]) / niterations
  write.csv(mala_iteration_time, sprintf("mala_time_d=%d.csv", dimension))
  
  mala_lag <- 10000
  nsamples <- 50
  max_iterations <- 100000
  meeting_times <- foreach(irep = 1:nsamples, .combine = rbind) %dorng% {
    res <- simulate_meeting_time(mala_single_kernel, mala_coupled_kernel, rinit, max_iterations = max_iterations, L = mala_lag)
    res$meeting_time
  }
  write.csv(meeting_times, sprintf("mala_meetings_d=%d_lag=%d.csv", dimension, mala_lag))
}

for (dimension in c(50, 100, 200, 300, 400, 500, 600, 800, 1000)) {
  run_for_dimension(dimension)
}
