# load packages
library(debiasedmcmc)
library(debiasedhmc)
library(doParallel)
library(doRNG)
#setmytheme()
rm(list = ls())
set.seed(21)
registerDoParallel(cores = detectCores()-2)
#


# target is Multivariate Normal
dimension <- 800
# define covariance matrix of target
Sigma_pi <- diag(1, dimension, dimension)
alpha <- 0.5
for (i in 1:dimension){
  for (j in 1:dimension){
    Sigma_pi[i,j] <- alpha^(abs(i-j))
  }
}

Sigma_pi_chol <- chol(Sigma_pi)
inv_Sigma_pi_chol <- solve(Sigma_pi_chol)
inv_Sigma_pi <- solve(Sigma_pi)

# target log-density
# target <- function(x) fast_dmvnorm(matrix(x, nrow = 1), rep(0, dimension), Sigma_pi)
target <- function(x) fast_dmvnorm_chol_inverse(matrix(x, nrow = 1), rep(0, dimension), inv_Sigma_pi_chol)

# gradient of target log-density
gradtarget <- function(x)- t(inv_Sigma_pi %*% x)[1,]


# proposal covariance matrix
Sigma_proposal <- diag(1, dimension, dimension)
Sigma_chol <- chol(Sigma_proposal)
inv_Sigma_chol <- solve(Sigma_chol)
#

# initial distribution
rinit <- function() fast_rmvnorm(1, rep(1, dimension), diag(dimension))[1,]

## Random walk Metropolis-Hastings with Normal increments

# Markov kernel of the chain
mh_single_kernel <- function(chain_state, current_pdf, stepsize){
  proposal_value <- fast_rmvnorm_chol(1, chain_state, stepsize * Sigma_chol)
  proposal_pdf <- target(proposal_value)
  accept <- (log(runif(1)) < (proposal_pdf - current_pdf))
  if (accept){
    return(list(chain_state = proposal_value, current_pdf = proposal_pdf))
  } else {
    return(list(chain_state = chain_state, current_pdf = current_pdf))
  }
}

#################################################################################### Markov kernel of the coupled chain
mh_coupled_kernel <- function(chain_state1, chain_state2, current_pdf1, current_pdf2, stepsize){
  proposal_value <- debiasedmcmc:::rnorm_reflectionmax_(chain_state1, chain_state2, stepsize * Sigma_chol, inv_Sigma_chol / stepsize)
  proposal1 <- proposal_value$xy[,1]
  proposal2 <- proposal_value$xy[,2]
  if (proposal_value$identical){
    proposal2 <- proposal1
  }
  proposal_pdf1 <- target(proposal1)
  proposal_pdf2 <- target(proposal2)
  logu <- log(runif(1))
  accept1 <- FALSE
  accept2 <- FALSE
  if (is.finite(proposal_pdf1)){
    accept1 <- (logu < (proposal_pdf1 - current_pdf1))
  }
  if (is.finite(proposal_pdf2)){
    accept2 <- (logu < (proposal_pdf2 - current_pdf2))
  }
  if (accept1){
    chain_state1 <- proposal1
    current_pdf1 <- proposal_pdf1
  }
  if (accept2){
    chain_state2 <- proposal2
    current_pdf2 <- proposal_pdf2
  }
  return(list(chain_state1 = chain_state1, chain_state2 = chain_state2, current_pdf1 = current_pdf1, current_pdf2 = current_pdf2))
}


# Unadjusted Langevin Algorithm (ULA)

ula_single_kernel <- function(chain_state, current_pdf, stepsize){
  mean_proposal <- chain_state + stepsize/2 * t(Sigma_chol %*%  gradtarget(chain_state))[1,]
  proposal_value <- fast_rmvnorm_chol(1, mean_proposal, sqrt(stepsize) * Sigma_chol)[1,]
  proposal_pdf <- target(proposal_value)
  return(list(chain_state = proposal_value, current_pdf = proposal_pdf))
}

ula_coupled_kernel <- function(chain_state1, chain_state2, current_pdf1, current_pdf2, stepsize){
  mean_proposal1 <- chain_state1 + stepsize/2 * t(Sigma_chol %*%  gradtarget(chain_state1))[1,]
  mean_proposal2 <- chain_state2 + stepsize/2 * t(Sigma_chol %*%  gradtarget(chain_state2))[1,]
  proposal_value <- debiasedmcmc:::rnorm_reflectionmax_(mean_proposal1, mean_proposal2, sqrt(stepsize) * Sigma_chol, inv_Sigma_chol / sqrt(stepsize))
  proposal1 <- proposal_value$xy[,1]
  proposal2 <- proposal_value$xy[,2]
  proposal_pdf1 <- target(proposal1)
  proposal_pdf2 <- target(proposal2)
  return(list(chain_state1 = proposal1, chain_state2 = proposal2, current_pdf1 = proposal_pdf1, current_pdf2 = proposal_pdf2))
}

## Metropolis adjusted Langevin Algorithm (MALA)
mala_single_kernel <- function(chain_state, current_pdf, stepsize){
  mean_proposal <- chain_state + stepsize/2 * t(Sigma_chol %*%  gradtarget(chain_state))[1,]
  proposal_value <- fast_rmvnorm_chol(1, mean_proposal, sqrt(stepsize) * Sigma_chol)[1,]
  mean_revert <- proposal_value + stepsize/2  * t(Sigma_chol %*% gradtarget(proposal_value))[1,]
  proposal_pdf <- target(proposal_value)
  accept <- (log(runif(1)) < (proposal_pdf - current_pdf + fast_dmvnorm_chol_inverse(matrix(chain_state, nrow = 1), mean_revert, inv_Sigma_chol / sqrt(stepsize)) -
                                fast_dmvnorm_chol_inverse(matrix(proposal_value, nrow = 1), mean_proposal, inv_Sigma_chol / sqrt(stepsize))))
  if (accept){
    return(list(chain_state = proposal_value, current_pdf = proposal_pdf))
  } else {
    return(list(chain_state = chain_state, current_pdf = current_pdf))
  }
}

mala_coupled_kernel <- function(chain_state1, chain_state2, current_pdf1, current_pdf2, stepsize){
  mean_proposal1 <- chain_state1 + stepsize/2 * t(Sigma_chol %*%  gradtarget(chain_state1))[1,]
  mean_proposal2 <- chain_state2 + stepsize/2 * t(Sigma_chol %*%  gradtarget(chain_state2))[1,]
  proposal_value <- debiasedmcmc:::rnorm_reflectionmax_(mean_proposal1, mean_proposal2, sqrt(stepsize) * Sigma_chol, inv_Sigma_chol / sqrt(stepsize))
  proposal1 <- proposal_value$xy[,1]
  proposal2 <- proposal_value$xy[,2]
  proposal_pdf1 <- target(proposal1)
  proposal_pdf2 <- target(proposal2)

  mean_revert1 <- proposal1 + stepsize/2  * t(Sigma_chol %*% gradtarget(proposal1))[1,]
  mean_revert2 <- proposal2 + stepsize/2  * t(Sigma_chol %*% gradtarget(proposal2))[1,]

  logu <- log(runif(1))

  accept1 <- FALSE
  accept2 <- FALSE
  if (is.finite(proposal_pdf1)){
    accept1 <- (logu < (proposal_pdf1 - current_pdf1 + fast_dmvnorm_chol_inverse(matrix(chain_state1, nrow = 1), mean_revert1, inv_Sigma_chol / sqrt(stepsize)) -
                          fast_dmvnorm_chol_inverse(matrix(proposal1, nrow = 1), mean_proposal1, inv_Sigma_chol / sqrt(stepsize))))
  }
  if (is.finite(proposal_pdf2)){
    accept2 <- (logu < (proposal_pdf2 - current_pdf2 + fast_dmvnorm_chol_inverse(matrix(chain_state2, nrow = 1), mean_revert2, inv_Sigma_chol / sqrt(stepsize)) -
                          fast_dmvnorm_chol_inverse(matrix(proposal2, nrow = 1), mean_proposal2, inv_Sigma_chol / sqrt(stepsize))))
  }
  if (accept1){
    chain_state1 <- proposal1
    current_pdf1 <- proposal_pdf1
  }
  if (accept2){
    chain_state2 <- proposal2
    current_pdf2 <- proposal_pdf2
  }
  return(list(chain_state1 = chain_state1, chain_state2 = chain_state2, current_pdf1 = current_pdf1, current_pdf2 = current_pdf2))
}

## Generic way of sampling meeting times with an arbitrary lag
rmeeting <- function(rinit, target, single_kernel, coupled_kernel, stepsize = 1, lag = 1, max_iterations = Inf){
  chain_state1 <- rinit()
  chain_state2 <- rinit()
  current_pdf1 <- target(chain_state1)
  current_pdf2 <- target(chain_state2)
  iter <- 0
  for (s in 1:lag){
    sres1 <- single_kernel(chain_state1, current_pdf1, stepsize)
    chain_state1 <- sres1$chain_state
    current_pdf1 <- sres1$current_pdf
    iter <- iter + 1
  }
  meet <- FALSE
  meetingtime <- Inf
  # iter here is 1; at this point we have X_1,Y_0 and we are going to generate successively X_t,Y_{t-1} where iter = t
  while (is.infinite(meetingtime) && iter < max_iterations){
    iter <- iter + 1
    res_coupled_kernel <- coupled_kernel(chain_state1, chain_state2, current_pdf1, current_pdf2, stepsize)
    chain_state1 <- res_coupled_kernel$chain_state1
    chain_state2 <- res_coupled_kernel$chain_state2
    current_pdf1 <- res_coupled_kernel$current_pdf1
    current_pdf2 <- res_coupled_kernel$current_pdf2
    if (all(chain_state1 == chain_state2) && !meet){
      # recording meeting time tau
      meet <- TRUE
      meetingtime <- iter
    }
  }
  return(list(meetingtime = meetingtime, iteration = iter))
}


### Test algorithms

#### RWMH
## Long run
niterations <- 1000
chain <- matrix(nrow = niterations, ncol = dimension)
chain_state <- rinit()
current_pdf <- target(chain_state)
ptm_irep <- proc.time()
for (iteration in 1:niterations){
  sres <- mh_single_kernel(chain_state, current_pdf, stepsize = 2/sqrt(dimension))
  chain_state <- sres$chain_state
  current_pdf <- sres$current_pdf
  chain[iteration,] <- chain_state
}
elapsedtime_rwmh <- as.numeric((proc.time() - ptm_irep)[3])/niterations
matplot(chain[1:1e4], type = "l")

hist(chain[1e3:niterations,1], prob = TRUE, nclass = 5e2)
curve(dnorm(x, mean = 0, sd = sqrt(Sigma_pi[1,1])), add = TRUE, lty = 2)
cov(chain[1e3:niterations,])
colMeans(chain[1e3:niterations,])





# Meeting times
nsamples <- 100
#meetings <- foreach(irep = 1:nsamples) %dorng% {
#  rmeeting(rinit, target, mh_single_kernel, mh_coupled_kernel, stepsize = 2/sqrt(dimension), lag = 1)
#}
#hist(sapply(meetings, function(x) x$meetingtime), main = "Meeting times with RWMH")

max_iterations <- 20000
t_start <- 0
t_end <- 1000
meetings_lag_15000 <- rep(0, nsamples)

for(irep in 1:nsamples){
  lag_15000 <- rmeeting(rinit, target, mh_single_kernel, mh_coupled_kernel, stepsize = 2/sqrt(dimension), lag = 15000)
  meetings_lag_15000[irep] <- lag_15000$meetingtime
  print(irep)
}

#write.csv(meetings_lag_15000, "rwmh_meetings_dim50_lag_15000.csv")

#meetingtimes_15000 <-  foreach(irep = 1:nsamples, .combine = c) %dorng% {
#  rmeeting(rinit, target, mh_single_kernel, mh_coupled_kernel, stepsize = 2/sqrt(dimension), lag = 15000)$meetingtime
#}



#### ULA
## Long run
niterations <- 100
chain <- matrix(nrow = niterations, ncol = dimension)
chain_state <- rinit()
current_pdf <- target(chain_state)
ptm_irep <- proc.time()
for (iteration in 1:niterations){
  sres <- ula_single_kernel(chain_state, current_pdf, 0.1 * (1/dimension)^(1/6))
  chain_state <- sres$chain_state
  current_pdf <- sres$current_pdf
  chain[iteration,] <- chain_state
}
elapsedtime_ula <- as.numeric((proc.time() - ptm_irep)[3])/niterations

# matplot(chain[1:1e3,2], type = "l")
# matplot(chain[1:1e4,2], type = "l")

hist(chain[1e3:niterations,1], prob = TRUE, nclass = 4e2)
curve(dnorm(x, mean = 0, sd = sqrt(Sigma_pi[1,1])), add = TRUE, lty = 2, col = "red", lwd = 2)
cov(chain[1e4:niterations,])
colMeans(chain[1e3:niterations,])

# ## Meeting times
nsamples <- 50
#meetings <- foreach(irep = 1:NREP) %dorng% {
#  rmeeting(rinit, target, ula_single_kernel, ula_coupled_kernel, stepsize = 0.1 * (1/dimension)^(1/6), lag = 1)
#}
#hist(sapply(meetings, function(x) x$meetingtime), main = "Meeting times with ULA")

max_iterations <- 100000
t_start <- 0
t_end <- 1000
meetings_lag_10000 <- rep(0, nsamples)

for(irep in 1:nsamples){
  lag_10000 <- rmeeting(rinit, target, ula_single_kernel, ula_coupled_kernel, stepsize = 0.1 * (1/dimension)^(1/6), lag = 10000)
  meetings_lag_10000[irep] <- lag_10000$meetingtime
  print(irep)
}

#write.csv(meetings_lag_10000, "ula_meetings_dim800_lag_10000.csv")


#### MALA
## Long run
niterations <- 100
chain <- matrix(nrow = niterations, ncol = dimension)
chain_state <- rinit()
current_pdf <- target(chain_state)
ptm_irep <- proc.time()
for (iteration in 1:niterations){
  sres <- mala_single_kernel(chain_state, current_pdf, (1/dimension)^(1/6))
  chain_state <- sres$chain_state
  current_pdf <- sres$current_pdf
  chain[iteration,] <- chain_state
}
elapsedtime_mala <- as.numeric((proc.time() - ptm_irep)[3])/niterations

matplot(chain[1:1e3,2], type = "l")
# matplot(chain[1:1e4,2], type = "l")
#
hist(chain[1e3:niterations,1], prob = TRUE, nclass = 4e2)
curve(dnorm(x, mean = 0, sd = sqrt(Sigma_pi[1,1])), add = TRUE, lty = 2, col = "red", lwd = 2)
cov(chain[1e4:niterations,])
colMeans(chain[1e3:niterations,])

## Meeting times
nsamples <- 50
#meetings <- foreach(irep = 1:NREP) %dorng% {
#  rmeeting(rinit, target, mala_single_kernel, mala_coupled_kernel, stepsize = (1/dimension)^(1/6), lag = 1)
#}
#hist(sapply(meetings, function(x) x$meetingtime), main = "Meeting times with MALA")

max_iterations <- 100000
t_start <- 0
t_end <- 5000
meetings_lag_10000 <- rep(0, nsamples)

for(irep in 1:nsamples){
  lag_10000 <- rmeeting(rinit, target, mala_single_kernel, mala_coupled_kernel, stepsize = (1/dimension)^(1/6), lag = 10000)
  meetings_lag_10000[irep] <- lag_10000$meetingtime
  print(irep)
}

#write.csv(meetings_lag_10000, "mala_meetings_dim800_lag_10000.csv")
  