library(unbiasedmcmc)
library(doParallel)
library(unbiasedmcmc)

source("coupling.R")

set.seed(1)
registerDoParallel(cores = detectCores())

# Load german credit dataset.
data(germancredit)
X <- scale(X)
X[,1] <- rep(1, nrow(X))
n <- nrow(X)
p <- ncol(X)
design_matrix <- unname(X)
tdesign_matrix <- t(design_matrix)
response <- Y
new_response <- 2*response - 1
nsamples <- nrow(design_matrix)
dimension <- ncol(design_matrix)

# Prior variance for regression coefficients.
sigma2 <- 10

stable_log_sigmoid <- function(x){
  output <- vector(mode = "logical", length = length(x))
  mask <- x > 0
  nmask <- !mask
  output[mask] <- -log(1+exp(-x[mask]))
  output[nmask] <- x[nmask] - log1p(exp(x[nmask]))
  return(output)
}

sigmoid <- function(x) {
  1 / (1 + exp(-x))
}

logtarget <- function(beta){
  xbeta <- design_matrix %*% beta
  loglikelihood <- sum(stable_log_sigmoid(new_response*xbeta))
  return(loglikelihood - sum(beta^2)/(2*sigma2))
}

gradlogtarget <- function(beta){
  xbeta <- design_matrix %*% beta
  tdesign_matrix %*% (sigmoid(-new_response * xbeta) * new_response) - beta / (2 * sigma2)
}

# Distribution of X_0.
rinit <- function() {
  x <- rnorm(dimension) * sqrt(sigma2)
  return(list(
    chain_state=x,
    current_pdf=logtarget(x)
  ))
}

nsamples <- 100
max_iterations <- 10000

hmc_experiments_df <- data.frame(
  stepsize = c(0.025, 0.025, 0.025, 0.025),
  nsteps = c(4, 5, 6, 7),
  lag = c(1000, 1000, 500, 2000)
)

for (i in 1:nrow(hmc_experiments_df)) {
  row <- hmc_experiments_df[i, ]
  filename = sprintf("hmc_meetings_lag=%d_stepsize=%f_nsteps=%d.csv", row$lag, row$stepsize, row$nsteps)
  print(filename)
  
  hmc <- get_hmc_kernel(logtarget, gradlogtarget, row$stepsize, row$nsteps, dimension)

  omega <- 1 / 20 # probability of selecting coupled RWMH
  Sigma_std <- 1e-3
  Sigma_proposal <- Sigma_std^2 * diag(1, dimension, dimension)
  mh <- get_mh_kernel(logtarget, Sigma_proposal, dimension)
  
  mixture_kernel <- function(chain_state, current_pdf, iteration) {
    if (runif(1) < omega){
      return(mh$kernel(chain_state, current_pdf, iteration))
    } else {
      return(hmc$kernel(chain_state, current_pdf, iteration))
    }
  }
  
  mixture_coupled_kernel <- function(chain_state1, chain_state2, current_pdf1, current_pdf2, iteration) {
    if (runif(1) < omega){
      return(mh$coupled_kernel(chain_state1, chain_state2, current_pdf1, current_pdf2, iteration))
    } else {
      return(hmc$coupled_kernel(chain_state1, chain_state2, current_pdf1, current_pdf2, iteration))
    }
  }
  
  meetings <- times(nsamples) %dopar% {
      res <- simulate_meeting_time(mixture_kernel, mixture_coupled_kernel, rinit, max_iterations = max_iterations, L = row$lag)
      res$meeting_time
    }
  
  write.csv(meetings, filename)
}
