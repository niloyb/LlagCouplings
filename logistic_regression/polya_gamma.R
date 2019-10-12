library(BayesLogit)
library(doParallel)
library(unbiasedmcmc)

source("coupling.R")

set.seed(1)
registerDoParallel(cores = detectCores())

## Loading data
# load german credit dataset
data(germancredit)
X <- scale(X)
X[,1] <- rep(1, nrow(X))
n <- nrow(X)
p <- ncol(X)

# Prior
b <- matrix(0, nrow = p, ncol = 1)
B <- diag(10, p, p)
logistic_setting <- logisticregression_precomputation(Y, X, b, B)

# Functions

rinit <- function(){
  x <- t(fast_rmvnorm(1, mean = b, covariance = B))
  return(list(
    chain_state=x,
    current_pdf=0
  )) 
}

single_kernel <- function(chain_state, current_pdf, iteration) {
  zs <- abs(logisticregression_xbeta(logistic_setting$X, t(chain_state)))
  w <- BayesLogit::rpg(logistic_setting$n, h=1, z=zs)
  res <- logisticregression_m_and_sigma(w, X, logistic_setting$invB, logistic_setting$KTkappaplusinvBtimesb)
  chain_state <- t(fast_rmvnorm_chol(1, res$m, res$Cholesky))
  return(list(
    chain_state=chain_state,
    current_pdf=0
  ))
}

coupled_kernel <- function(chain_state1, chain_state2, current_pdf1, current_pdf2, iteration) {
  ws <- sample_w(chain_state1, chain_state2, logistic_setting$X, mode='rej_samp')
  betas <- sample_beta(ws$w1, ws$w2, logistic_setting)
  if (all(ws$w1 == ws$w2)){
    betas$beta2 <- betas$beta1
  }
  return(list(
    chain_state1=cbind(betas$beta1),
    chain_state2=cbind(betas$beta2),
    current_pdf1=0,
    current_pdf2=0
  ))
}

# Generating meeting times
#nsamples <- 1000
nsamples <- 100
max_iterations <- 20000
lag <- 350
meetings <- times(nsamples) %dopar% {
  res <- simulate_meeting_time(single_kernel, coupled_kernel, rinit, max_iterations, L = lag)
  res$meeting_time
}

write.csv(meetings, sprintf("PG_meetings_lag=%d.csv", lag))
