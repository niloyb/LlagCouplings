
# devtools::install_github('pierrejacob/debiasedhmc')
rm(list=ls())
library(tictoc)
library(debiasedmcmc)
library(debiasedhmc)
library(coda)

source("tv_wasserstein_bounds.R")

# Loading data
#load german credit dataset
data(germancredit)
X <- scale(X)
X[,1] <- rep(1, nrow(X))
n <- nrow(X)
p <- ncol(X)
design_matrix <- X
response <- Y
new_response <- 2*response - 1
nsamples <- nrow(design_matrix)
dimension <- ncol(design_matrix)


# Model for HMC
sigma2 <- 10

stable_log_sigmoid <- function(x){
  output <- -log(1+exp(-x))
  output[x<0] <- x[x<0] - log1p(exp(x[x<0]))
  return(output)
}

logtarget <- function(beta){
  xbeta <- crossprod(t(design_matrix), beta)
  loglikelihood <- sum(stable_log_sigmoid(new_response*xbeta))
  return(loglikelihood - sum(beta^2)/(2*sigma2))
}

gradlogtarget <- function(beta){
  xbeta <- crossprod(t(design_matrix), beta)
  grad_log <- (crossprod((design_matrix), pracma::sigmoid(-new_response*xbeta)*new_response) - beta/(2*sigma2))[,1]
  grad_log <- unname(grad_log)
  return( grad_log )
}


# initial distribution
rinit <- function() rnorm(dimension)

################## Parameters 
# optimal HMC parameters
# stepsize <- 0.000025
stepsize <- 0.025
nsteps <- 6

################## HMC kernels and coupled kernels
# define HMC kernel and coupled HMC kernel
hmc <- get_hmc_kernel(logtarget, gradlogtarget, stepsize, nsteps, dimension)

# define RWMH kernel and coupled RWMH kernel
omega <- 1 / 20 # probability of selecting coupled RWMH
Sigma_std <- 1e-3 # proposal standard deviation of RWMH
Sigma_proposal <- Sigma_std^2 * diag(1, dimension, dimension)
mh <- get_mh_kernel(logtarget, Sigma_proposal, dimension)

# define mixture kernel
mixture_kernel <- function(chain_state, current_pdf, iteration){
  if (runif(1) < omega){
    return(mh$kernel(chain_state, current_pdf, iteration))
  } else {
    return(hmc$kernel(chain_state, current_pdf, iteration))
  }
}

# define coupled mixture kernel
mixture_coupled_kernel <- function(chain_state1, chain_state2, current_pdf1, current_pdf2, iteration){
  if (runif(1) < omega){
    return(mh$coupled_kernel(chain_state1, chain_state2, current_pdf1, current_pdf2, iteration))
  } else {
    return(hmc$coupled_kernel(chain_state1, chain_state2, current_pdf1, current_pdf2, iteration))
  }
}


# Behavior from a single run
# current_value1 <- rinit()
# pdf_value1 <- logtarget(current_value1)
# current_value2 <- rinit()
# pdf_value2 <- logtarget(current_value2)
# 
# niterations <- 500
# chain1_beta <- chain2_beta <- matrix(ncol=p, nrow=niterations)
# chain1_beta[1,] <- current_value1
# chain2_beta[1,] <- current_value2
# 
# ptc <- proc.time()
# for(t in 2:niterations){
#   current_value <- mixture_coupled_kernel(current_value1, current_value2, pdf_value1, pdf_value2)
#   chain1_beta[t,] <- current_value1 <- current_value$chain_state1
#   chain2_beta[t,] <- current_value2 <- current_value$chain_state2
#   pdf_value1 <- current_value$current_pdf1
#   pdf_value2 <- current_value$current_pdf2
#   print(t)
#   #if(all(current_value1==current_value2)) break
#   # current_value <- hmc$kernel(current_value1, pdf_value1)
#   # print(current_value$chain_state[1])
#   # chain1_beta[t,] <- current_value1 <- current_value$chain_state
#   # pdf_value1 <- current_value$current_pdf
#   # print(current_value$accept)
# }
# meetingtime=t
# time_taken <- (proc.time()[3]-ptc[3])/t
# 
# matplot(chain1_beta, type = 'l')
# matplot(chain2_beta, type = 'l')
# 
# plot(sapply(1:niterations, function(index) sum((chain1_beta[index,] - chain2_beta[index,])^2)), type = "l", log = c("y"))



################################################################## Meeting times HMC
mh_hmc_meeting <- function(single_kernel, coupled_kernel, rinit, max_iterations = Inf, L=1, save_chains=FALSE){
  
  chain1 <- NA
  chain2 <- NA
  
  if (save_chains==TRUE)
  {
    if(max_iterations != Inf)
    {
      chain1 <- matrix(0, nrow = max(max_iterations,L), ncol = p)
      chain2 <- matrix(0, nrow = max(max_iterations,L), ncol = p)
    } else {
      chain1 <- matrix(0, nrow = 1000, ncol = p)
      chain2 <- matrix(0, nrow = 1000, ncol = p)
    }
  }
  
  chain_state1 <- rinit()
  current_pdf1 <- logtarget(chain_state1)
  chain_state2 <- rinit()
  current_pdf2 <- logtarget(chain_state2)
  
  if (L>1)
  {
    for (l in 1:(L-1))
    {
      output <- single_kernel(chain_state1, current_pdf1)
      chain_state1 <- output$chain_state
      current_pdf1 <- output$current_pdf
    }
  }
  
  if (save_chains==TRUE){
    chain1[1,] <- chain_state1
    chain2[1,] <- chain_state2
  }
  
  output <- single_kernel(chain_state1, current_pdf1)
  chain_state1 <- output$chain_state
  current_pdf1 <- output$current_pdf
  
  if (save_chains==TRUE){
    chain1[2,] <- chain_state1
  }
  
  iter <- L
  meet <- FALSE
  finished <- FALSE
  meetingtime <- Inf
  # iter here is L; at this point we have X_L,Y_0 and we are going to generate successively X_t,Y_{t-1} where iter = t
  while (!finished && iter < max_iterations){
    iter <- iter + 1
    if (meet){
      output <- single_kernel(chain_state1, current_pdf)
      chain_state1 <- output$chain_state
      current_pdf1 <- output$current_pdf
      chain_state2 <- chain_state1
      current_pdf2 <- current_pdf1
    } else {
      res_coupled_kernel <- coupled_kernel(chain_state1, chain_state2, current_pdf1, current_pdf2)
      chain_state1 <- res_coupled_kernel$chain_state1
      current_pdf1 <- res_coupled_kernel$current_pdf1
      chain_state2 <- res_coupled_kernel$chain_state2
      current_pdf2 <- res_coupled_kernel$current_pdf2
      if (all(chain_state1 == chain_state2) && !meet){
        # recording meeting time tau
        meet <- TRUE
        meetingtime <- iter
      }
    }
    
    
    
    if (save_chains==TRUE) {
      chain1[((iter)+2-L),] <- chain_state1
      chain2[((iter)+1-L),] <- chain_state2
      
      # Appending extra memory for chains
      if ((iter)+2-L+1000> dim(chain1)[1] )
      {
        chain1 <- rbind(chain1, matrix(0, nrow = 1000, ncol = p))
        chain2 <- rbind(chain2, matrix(0, nrow = 1000, ncol = p))
      }
      
    }
    
    # stop after max(m, tau) steps
    if (iter >= meetingtime){
      finished <- TRUE
    }
    
    
  }
  return(list(meetingtime = meetingtime, chain1 = chain1, chain2 = chain2, iteration = iter, finished = finished))
}


nsamples <- 100
max_iterations <- 10000
t_start <- 0
t_end <- 1000
meetings_lag_2000 <- rep(0, nsamples)

for(irep in 1:nsamples){
  lag_2000 <- mh_hmc_meeting(mixture_kernel, mixture_coupled_kernel, rinit, max_iterations, L = 2000)
  meetings_lag_2000[irep] <- lag_2000$meetingtime
  
  print(irep)
}


#write.csv(meetings_lag_2000, "hmc_meetings_lag_2000_stepsize_0.025_nsteps6.csv")


