rm(list = ls())

library(doParallel)
library(doRNG)

# generate meeting time
simulate_meeting_time <-
  function(single_kernel,
           coupled_kernel,
           rinit,
           max_iterations = Inf,
           L = 1) {
    # initialize first chain
    init_res1 <- rinit()
    chain_state1 <- init_res1$chain_state
    current_pdf1 <- init_res1$current_pdf
    # advance for L steps to create a lag
    for (t in 1:L) {
      output <- single_kernel(chain_state1, current_pdf1)
      chain_state1 <- output$chain_state
      current_pdf1 <- output$current_pdf
    }
    # initialize second chain
    init_res2 <- rinit()
    chain_state2 <- init_res2$chain_state
    current_pdf2 <- init_res2$current_pdf
    # then propagate coupled chains until they meet
    t <- L + 1
    meeting_time <- Inf
    while (t <= max_iterations) {
      res_coupled_kernel <-
        coupled_kernel(chain_state1, chain_state2, current_pdf1, current_pdf2)
      chain_state1 <- res_coupled_kernel$chain_state1
      current_pdf1 <- res_coupled_kernel$current_pdf1
      chain_state2 <- res_coupled_kernel$chain_state2
      current_pdf2 <- res_coupled_kernel$current_pdf2
      t <- t + 1
      if (res_coupled_kernel$identical) {
        # record meeting time
        meeting_time <- t - 1
        break
      }
    }
    return(list(
      meeting_time = meeting_time
    ))
  }

## compute TV upper bound at time t 
## from meeting times generated with lag L
tv_upper_bound_estimates <- function(meeting_times, L, t){
  return(pmax(0,ceiling((meeting_times-L-t)/L)))
}

## import source files from the github repository
source("https://raw.githubusercontent.com/niloyb/LlagCouplings/master/ula_versus_mala/helper_functions/mvnorm.R")
source("https://raw.githubusercontent.com/niloyb/LlagCouplings/master/ula_versus_mala/helper_functions/mvnorm_couplings.R")
## import c++ file and compile 
download.file(url="https://raw.githubusercontent.com/niloyb/LlagCouplings/master/ula_versus_mala/helper_functions/mvnorm.h", destfile="mvnorm.h")
mvnormcppfile <- "https://raw.githubusercontent.com/niloyb/LlagCouplings/master/ula_versus_mala/helper_functions/mvnorm.cpp"
locfile <- "mvnorm.cpp"
download.file(url=mvnormcppfile, destfile=locfile)
Rcpp::sourceCpp(locfile)

## set random seed
set.seed(1)
## register cores
registerDoParallel(cores = detectCores()-1)

## function that generates meeting times for ULA, MALA
## nsamples = number of meeting times, generated in parallel
## step_size_cnst = constant in front of n^{-1/6} that defines the stepsize 
## ula_lag and mala_lag refer to the lag employed to generate meeting times

run_ULAMALA <- function(dimension, nsamples=50, step_size_cnst=1,
                        ula_lag = 1e3, mala_lag = 1e3) {
  # define covariance matrix of Normal target
  Sigma_pi <- diag(1, dimension, dimension)
  alpha <- 0.5
  for (i in 1:dimension) {
    for (j in 1:dimension) {
      Sigma_pi[i, j] <- alpha^(abs(i - j))
    }
  }
  # precomputation of Cholesky factor of Sigma and its inverse
  Sigma_pi_chol <- chol(Sigma_pi)
  inv_Sigma_pi_chol <- solve(Sigma_pi_chol)
  inv_Sigma_pi <- solve(Sigma_pi)
  # target log-density
  target <- function(x) {
    fast_dmvnorm_chol_inverse(matrix(x, nrow = 1), rep(0, dimension), inv_Sigma_pi_chol)
  }
  # gradient of target log-density
  gradtarget <- function(x) {
    -(inv_Sigma_pi %*% x)[,1]
  }
  # proposal covariance matrix (to be multiplied by the squared stepsize)
  Sigma_proposal <- diag(1, dimension, dimension)
  # precomputation of Cholesky factor of Sigma_proposal and its inverse
  Sigma_chol <- chol(Sigma_proposal)
  inv_Sigma_chol <- solve(Sigma_chol)
  # initial distribution Normal(1, I)
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
  # stepsize as a function of dimension
  ula_stepsize <- step_size_cnst * (dimension^(-1/6))
  # ULA transition kernel
  ula_single_kernel <- function(chain_state, current_pdf) {
    mean_proposal <- chain_state + (ula_stepsize)^2 / 2 * t(Sigma_proposal %*% gradtarget(chain_state))[1,]
    # note the function 'fast_rmvnorm_chol' is parameterized with the Cholesky factor of the variance
    proposal_value <- fast_rmvnorm_chol(1, mean_proposal, ula_stepsize * Sigma_chol)[1,]
    return(list(
      chain_state = proposal_value,
      current_pdf = NA
    ))
  }
  # coupled ULA transition kernel
  ula_coupled_kernel <- function(chain_state1, chain_state2, current_pdf1, current_pdf2) {
    mean_proposal1 <- chain_state1 + (ula_stepsize)^2/2 * t(Sigma_proposal %*% gradtarget(chain_state1))[1,]
    mean_proposal2 <- chain_state2 + (ula_stepsize)^2/2 * t(Sigma_proposal %*% gradtarget(chain_state2))[1,]
    # note the function 'rmvnorm_reflectionmax' is parameterized with the Cholesky factor of the variance
    proposal_value <- rmvnorm_reflectionmax(mean_proposal1, mean_proposal2, ula_stepsize * Sigma_chol, inv_Sigma_chol / ula_stepsize)
    proposal1 <- proposal_value$xy[,1]
    proposal2 <- proposal_value$xy[,2]
    return(list(
      chain_state1 = proposal1,
      chain_state2 = proposal2,
      current_pdf1 = NA,
      current_pdf2 = NA,
      identical = proposal_value$identical
    ))
  }
  # time a single ULA chain
  niterations <- 10000
  res <- rinit()
  chain_state <- res$chain_state
  current_pdf <- target(chain_state)
  start_time <- proc.time()
  for (iteration in 1:niterations) {
    res <- ula_single_kernel(chain_state, current_pdf)
    chain_state <- res$chain_state
    current_pdf <- res$current_pdf
  }
  ula_iteration_time <- as.numeric((proc.time() - start_time)[3]) / niterations
  # simulate ULA meeting times
  ula_meetings <- foreach(irep = 1:nsamples, .combine = rbind) %dorng% {
    res <- simulate_meeting_time(ula_single_kernel, ula_coupled_kernel, rinit, L = ula_lag)
    res$meeting_time
  }
  ula_df <- data.frame(algo='ula', meetings=as.vector(ula_meetings), dimension=dimension, 
                       lag=ula_lag, iteration_time=ula_iteration_time, 
                       cnst = step_size_cnst,
                       stepsize=ula_stepsize, acceptance_prob=1)
  
  ###############################################
  # Metropolis adjusted Langevin Algorithm (MALA)
  ###############################################
  # stepsize as a function of dimension
  mala_stepsize = step_size_cnst * (dimension^(-1/6))
  # MALA transition kernel
  mala_single_kernel <- function(chain_state, current_pdf) {
    mean_proposal <- chain_state + (mala_stepsize)^2/2 * t(Sigma_proposal %*%  gradtarget(chain_state))[1,]
    # again check that 'fast_rmvnorm_chol' takes the sqrt of a variance as argument 
    proposal_value <- fast_rmvnorm_chol(1, mean_proposal, mala_stepsize * Sigma_chol)[1,]
    mean_revert <- proposal_value + (mala_stepsize)^2/2  * t(Sigma_proposal %*% gradtarget(proposal_value))[1,]
    proposal_pdf <- target(proposal_value)
    # fast_dmvnorm_chol_inverse is parametrized by the Cholesky factor of the inverse of proposal variance 
    accept <- (log(runif(1)) < (proposal_pdf - current_pdf + fast_dmvnorm_chol_inverse(matrix(chain_state, nrow = 1), mean_revert, inv_Sigma_chol / mala_stepsize) -
                                  fast_dmvnorm_chol_inverse(matrix(proposal_value, nrow = 1), mean_proposal, inv_Sigma_chol / mala_stepsize)))
    if (accept){
      return(list(chain_state = proposal_value, current_pdf = proposal_pdf, accept = accept))
    } else {
      return(list(chain_state = chain_state, current_pdf = current_pdf, accept = accept))
    }
  }
  # coupled MALA transition kernel
  mala_coupled_kernel <- function(chain_state1, chain_state2, current_pdf1, current_pdf2) {
    mean_proposal1 <- chain_state1 + (mala_stepsize)^2/2 * t(Sigma_proposal %*% gradtarget(chain_state1))[1,]
    mean_proposal2 <- chain_state2 + (mala_stepsize)^2/2 * t(Sigma_proposal %*% gradtarget(chain_state2))[1,]
    # note the function 'rmvnorm_reflectionmax' is parameterized with the Cholesky factor of the variance
    proposal_value <- rmvnorm_reflectionmax(mean_proposal1, mean_proposal2, mala_stepsize * Sigma_chol, inv_Sigma_chol / mala_stepsize)
    proposal1 <- proposal_value$xy[,1]
    proposal2 <- proposal_value$xy[,2]
    proposal_pdf1 <- target(proposal1)
    proposal_pdf2 <- target(proposal2)
    mean_revert1 <- proposal1 + (mala_stepsize)^2/2 * t(Sigma_proposal %*% gradtarget(proposal1))[1,]
    mean_revert2 <- proposal2 + (mala_stepsize)^2/2 * t(Sigma_proposal %*% gradtarget(proposal2))[1,]
    logu <- log(runif(1))
    accept1 <- FALSE
    accept2 <- FALSE
    if (is.finite(proposal_pdf1)){
      accept1 <- (logu < (proposal_pdf1 - current_pdf1 + fast_dmvnorm_chol_inverse(matrix(chain_state1, nrow = 1), mean_revert1, inv_Sigma_chol / mala_stepsize) -
                            fast_dmvnorm_chol_inverse(matrix(proposal1, nrow = 1), mean_proposal1, inv_Sigma_chol / mala_stepsize)))
    }
    if (is.finite(proposal_pdf2)){
      accept2 <- (logu < (proposal_pdf2 - current_pdf2 + fast_dmvnorm_chol_inverse(matrix(chain_state2, nrow = 1), mean_revert2, inv_Sigma_chol / mala_stepsize) -
                            fast_dmvnorm_chol_inverse(matrix(proposal2, nrow = 1), mean_proposal2, inv_Sigma_chol / mala_stepsize)))
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
      current_pdf2 = current_pdf2,
      identical = (proposal_value$identical && accept1 && accept2)
    ))
  }
  # time MALA chain and record acceptance rate
  niterations <- 10000
  res <- rinit()
  chain_state <- res$chain_state
  current_pdf <- target(chain_state)
  start_time <- proc.time()
  acceptances <- rep(NA, niterations)
  for (iteration in 1:niterations) {
    res <- mala_single_kernel(chain_state, current_pdf)
    acceptances[iteration] <- res$accept
    chain_state <- res$chain_state
    current_pdf <- res$current_pdf
  }
  mala_iteration_time <- as.numeric((proc.time() - start_time)[3]) / niterations
  acceptance_prob <- mean(acceptances)
  # generate MALA meeting times
  mala_meetings <- foreach(irep = 1:nsamples, .combine = rbind) %dorng% {
    res <- simulate_meeting_time(mala_single_kernel, mala_coupled_kernel, rinit, L = mala_lag)
    res$meeting_time
  }
  mala_df <- data.frame(algo='mala', meetings=as.vector(mala_meetings), dimension=dimension, 
                        lag=mala_lag, iteration_time=mala_iteration_time,
                        cnst = step_size_cnst,
                        stepsize=mala_stepsize, acceptance_prob=acceptance_prob)
  return(rbind(ula_df, mala_df))
}

######################## running simulations ########################
# run_ULAMALA(dimension = 10, nsamples = 50, step_size_cnst = 1)

## constants in front of stepsize
## warning: if the constant is too large then ULA might not converge
cnst_seq <- c(0.75, 1, 1.25, 1.5)
## dimensions to consider
dimension_seq <- c(50, 100, 150, 200, 250)
## number of meeting times
nsamples <- 50
## run everything
output <- data.frame()
for (cnst in cnst_seq){
  cat("cnst =", cnst, "\n")
  for (dimension in dimension_seq){
    cat("dimension =", dimension, "\n")
    output <- rbind(output, 
                    run_ULAMALA(dimension = dimension, nsamples = nsamples, step_size_cnst = cnst))
  }
}
head(output)

################################ Plots ################################
library(dplyr)
library(ggplot2)
library(latex2exp)

# acceptance rate of MALA as a function of stepsize constant & dimension

ardf <- output %>% filter(algo == 'mala') %>%
  select(dimension, cnst, acceptance_prob) %>% group_by(dimension, cnst) %>%
  summarise(ar = mean(acceptance_prob))

g_accept <- ggplot(ardf, aes(x = dimension, y = ar, linetype = factor(cnst))) + geom_line() +
  geom_point() + 
  ggtitle("Proposals with stepsize h = C*d^{-1/6}") +
  xlab("dimension") + 
  ylab(TeX("acceptance rate")) + 
  theme_classic(base_size = 16) +
  theme(legend.position = "bottom", legend.box = "horizontal") +
  labs(linetype = "C") + 
  ylim(0,1) + scale_x_continuous(breaks = dimension_seq)
g_accept

# ggsave(filename = "mala_acceptrate2.pdf", plot = g_accept,
#        width = 8, height = 4)

## upper bounds on mixing times
## i.e. iteration such that TV < delta
delta <- 0.25

## retrieve meeting times for each configuraton
## and compute upper bound on mixing time from meeting times
## also compute an adjusted version, taking wall time per iteration into account
mixingdf <- data.frame()
for (algo_ in c("ula", "mala")){ 
  for (cnst_ in cnst_seq){
    for (dimension_ in dimension_seq){
      ## extract meeting times
      meeting_times <- output %>% filter(algo == algo_, cnst == cnst_, dimension == dimension_) %>%
        pull(meetings)
      ## retrieve lag value
      lag_ <- (output %>% filter(algo == algo_, cnst == cnst_, dimension == dimension_) %>%
        pull(lag))[1]
      ## retrieve wall time per iteration
      iteration_time_ <- (output %>% filter(algo == algo_, cnst == cnst_, dimension == dimension_) %>%
                            pull(iteration_time))[1]
      ## construct upper bound on TV at time t for t=t_start,...,t_end
      t_start <- 0
      t_end <- 1000
      tv_ub_values <- sapply(t_start:t_end, function(t) mean(tv_upper_bound_estimates(meeting_times, L=lag_, t)))
      ## find first t such that upper bound at t is less than delta
      mixing_time <- which(tv_ub_values < delta)[1]
      ## multiply by wall time
      adjusted_mixing_time <- mixing_time * iteration_time_
      ## record results in data frame
      mixingdf <- rbind(mixingdf,
                        data.frame(algo = algo_, cnst = cnst_,
                                   dimension = dimension_,
                                   mixing_time = mixing_time,
                                   adjusted_mixing_time = adjusted_mixing_time))
    }
  }
}

head(mixingdf)

g_mixing <- ggplot(mixingdf, aes(x = dimension, y = mixing_time,
                     linetype = factor(algo),
                     colour = factor(cnst))) + geom_line() + geom_point()
g_mixing <- g_mixing + 
  xlab("dimension") + 
  ylab(TeX("$t_{mix}(0.25)$")) + 
  theme_classic(base_size = 16) +
  theme(legend.position = "bottom", legend.box = "horizontal") +
  labs(linetype = "algo", colour = "C") +
  ggtitle("Proposals with stepsize h = C*d^{-1/6}") +
  scale_x_continuous(breaks = dimension_seq) + guides(colour = "none")
g_mixing <- g_mixing + facet_wrap(~cnst, labeller = label_both, nrow = 2)

g_mixing 

# ggsave(filename = "ula_mala_mixing2.pdf", plot = g_mixing,
#        width = 8, height = 6)

