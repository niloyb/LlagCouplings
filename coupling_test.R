source("coupling.R")

rinit <- function() {
  return(list(
    state=1,
    logpdf=0
  ))
}

single_kernel <- function(state, logpdf, iteration) {
  return(list(
    state=2,
    logpdf=0
  ))
}

coupled_calls <- 0
coupled_kernel <- function(state1, state2, logpdf1, logpdf2, iteration) {
  t <- get("coupled_calls", envir=.GlobalEnv)
  t <- t + 1
  assign("coupled_calls", t, envir=.GlobalEnv)
  if (t > 5) {
      return(list(
        state1=5,
        state2=5
      ))
  } else {
    return(list(
      state1=3,
      state2=9,
      logpdf1=0,
      logpdf2=0
    ))
  }
}

coupled_calls <- 0
simulate_meeting_time(single_kernel, coupled_kernel, rinit, max_iterations=100, L=1)

coupled_calls <- 0
simulate_meeting_time(single_kernel, coupled_kernel, rinit, max_iterations=5, L=1)

coupled_calls <- 0
simulate_meeting_time(single_kernel, coupled_kernel, rinit, max_iterations=100, L=5)

coupled_calls <- 0
simulate_ipm_bounds(single_kernel, coupled_kernel, rinit, tv_ipm_bound, max_iterations = 1000, L = 1)
coupled_calls <- 0
simulate_ipm_bounds(single_kernel, coupled_kernel, rinit, tv_ipm_bound, max_iterations = 1000, L = 5)
coupled_calls <- 0
simulate_ipm_bounds(single_kernel, coupled_kernel, rinit, tv_ipm_bound, max_iterations = 1000, L = 10)

coupled_calls <- 0
simulate_ipm_bounds(single_kernel, coupled_kernel, rinit, wasserstein_ipm_bound, max_iterations = 1000, L = 1)
coupled_calls <- 0
simulate_ipm_bounds(single_kernel, coupled_kernel, rinit, wasserstein_ipm_bound, max_iterations = 1000, L = 5)
coupled_calls <- 0
simulate_ipm_bounds(single_kernel, coupled_kernel, rinit, wasserstein_ipm_bound, max_iterations = 1000, L = 10)

