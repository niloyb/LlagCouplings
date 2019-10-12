default_cmp <- function(chain_state1, chain_state2) {
  all(chain_state1 == chain_state2)
}

# Simulate a single meeting time.
simulate_meeting_time <-
  function(single_kernel,
           coupled_kernel,
           rinit,
           max_iterations = Inf,
           L = 1,
           cmp = default_cmp) {
    init_res1 <- rinit()
    chain_state1 <- init_res1$chain_state
    current_pdf1 <- init_res1$current_pdf
    
    for (t in 1:L) {
      output <- single_kernel(chain_state1, current_pdf1, t)
      chain_state1 <- output$chain_state
      current_pdf1 <- output$current_pdf
    }
    
    init_res2 <- rinit()
    chain_state2 <- init_res2$chain_state
    current_pdf2 <- init_res2$current_pdf
    
    t <- L + 1
    meeting_time <- Inf
    while (t <= max_iterations) {
      res_coupled_kernel <-
        coupled_kernel(chain_state1, chain_state2, current_pdf1, current_pdf2, t)
      chain_state1 <- res_coupled_kernel$chain_state1
      current_pdf1 <- res_coupled_kernel$current_pdf1
      chain_state2 <- res_coupled_kernel$chain_state2
      current_pdf2 <- res_coupled_kernel$current_pdf2

      t <- t + 1
      if (cmp(chain_state1, chain_state2)) {
        # recording meeting time tau
        meeting_time <- t - 1
        break
      }
    }
    return(list(
      meeting_time = meeting_time
    ))
  }

# Simulate an IPM bound.
simulate_ipm_bounds <-
  function(single_kernel,
           coupled_kernel,
           rinit,
           ipm_bound_fn,
           max_iterations = Inf,
           L = 1) {
    init_res1 <- rinit()
    chain_state1 <- init_res1$chain_state
    current_pdf1 <- init_res1$current_pdf
    
    for (t in 1:L) {
      output <- single_kernel(chain_state1, current_pdf1, t)
      chain_state1 <- output$chain_state
      current_pdf1 <- output$current_pdf
    }
    
    init_res2 <- rinit()
    chain_state2 <- init_res2$chain_state
    current_pdf2 <- init_res2$current_pdf
    
    bounds <- c(ipm_bound_fn(chain_state1, chain_state2))
    
    t <- L + 1
    meeting_time <- Inf
    while (t <= max_iterations) {
      res_coupled_kernel <-
        coupled_kernel(chain_state1, chain_state2, current_pdf1, current_pdf2, t)
      chain_state1 <- res_coupled_kernel$chain_state1
      current_pdf1 <- res_coupled_kernel$current_pdf1
      chain_state2 <- res_coupled_kernel$chain_state2
      current_pdf2 <- res_coupled_kernel$current_pdf2
      
      # Add the bound for X_t, Y_{t-L} to all t - jL, j >= 0.
      ipm_bound <- ipm_bound_fn(chain_state1, chain_state2)
      bounds <- c(bounds, 0)
      prev_t <- t - L
      while (prev_t >= 0) {
        bounds[prev_t + 1] <- bounds[prev_t + 1] + ipm_bound
        prev_t <- prev_t - L
      }
      
      if (all(chain_state1 == chain_state2)) {
        meeting_time <- t
        break
      }
      t <- t + 1
    }
    return(list(bounds = bounds,
                meeting_time = meeting_time))
  }

tv_ipm_bound <- function(chain_state1, chain_state2) {
  if (all(chain_state1 == chain_state2)) {
    return(0)
  } else {
    return(1)
  }
}

wasserstein_ipm_bound <- function(chain_state1, chain_state2) {
  return(sum(abs(chain_state1 - chain_state2)))
}
