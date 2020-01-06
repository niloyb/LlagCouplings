# This script reproduces Figure 1 in the article 
# "Estimating Convergence of Markov chains with L-Lag Couplings", 
# by Niloy Biswas, Pierre E. Jacob and Paul Vanetti. 
# In Advances in Neural Information Processing Systems 32 (2019), 7389--7399.
# https://papers.nips.cc/paper/8958-estimating-convergence-of-markov-chains-with-l-lag-couplings

library(doParallel)
library(doRNG)
library(ggplot2)
library(ggridges)
library(gridExtra)
library(latex2exp)
library(reshape2)
library(tidyverse)

source("coupling.R")
source("tv_wasserstein_bounds.R")

set.seed(1)
registerDoParallel(cores = detectCores())

# Simple random walk MH targetting N(0,1) and starting with X_0 = initial_value. 
# Proposals will be simulated from a N(x, sd_proposal).
initial_value <- 10
sd_proposal <- 0.5

################
# Exact results.
################
# For "exact" results, we will simulate many MCMC chains and 
# approximate the densities using those samples.
chain_length <- 160
number_of_chains <- 50000

all_chains <- foreach(i = 1:number_of_chains, .combine = rbind) %dorng% {
  X <- matrix(0, nrow = 1, ncol = chain_length)
  X[1] <- initial_value
  for (t in 2:chain_length) {
    proposal <- rnorm(1, mean = X[t - 1], sd = sd_proposal)
    if (log(runif(1)) < -0.5 * (proposal ^ 2 - X[t - 1] ^ 2)) {
      X[t] <- proposal
    } else {
      X[t] <- X[t - 1]
    }
  }
  X
}

# Exact TV calculation from kernel density estimates.
nmcmc <- chain_length
X <- t(all_chains)
target <- function(x) {
  evals <- dnorm(x, mean = 0, sd = 1, log = TRUE)
  return(max(evals) + log(sum(exp(evals - max(evals)))))
}

## If some numerical error arises, then the TV is set to zero
exact_tvs <- foreach(imcmc = 1:chain_length, .combine = rbind) %dopar% {
  # compute kernel density estimate, based on Markov chains
  kde_ <- density(X[imcmc, ])
  # compute |f(x) - pi(x)| where f is the kernel density estimate based on samples, and pi is the target
  diff_12 <- function(x) {
    sapply(x, function(v) {
      if ((v > min(kde_$x)) && (v < max(kde_$x))) {
        return(abs(approx(kde_$x, kde_$y, xout = v)$y - exp(target(v))))
      } else {
        return(abs(exp(target(v))))
      }
    })
  }
  # plot curve of |f(x) - pi(x)| against x
  # curve(diff_12(x), from = -5, to = 15, log = "y", main = paste("iteration", imcmc))
  # try to perform numerical integration
  result_ <-
    try(integrate(diff_12, lower = -5, upper = 15, stop.on.error = FALSE)$value / 2)
  if (inherits(result_, "try-error")) {
    result_ <- 0
  }
  result_
}

# Exact Wasserstein calculations.
exact_wasserstein <- rep(0, chain_length)
for (t in 1:chain_length) {
  std_norm_samples <- rnorm(number_of_chains)
  std_norm_samples <- sort(std_norm_samples)
  marginal_samples <- all_chains[, t]
  marginal_samples <- sort(marginal_samples)
  exact_wasserstein[t] <-
    mean(abs(std_norm_samples - marginal_samples))
}

#####################
# Coupling algorithm.
#####################
logtarget <- function(x) {
  dnorm(x, mean = 0, sd = 1, log = TRUE)
}

sample_initial <- function() {
  return(list(
    chain_state = initial_value,
    current_pdf = logtarget(initial_value)
  ))
}

singleMH_kernel <- function(state, logpdf, iteration) {
  proposal <- rnorm(1, mean = state, sd = sd_proposal)
  proposal_logpdf <- logtarget(proposal)
  logu <- log(runif(1))
  if (logu < (proposal_logpdf - logpdf)) {
    state <- proposal
    logpdf <- proposal_logpdf
  }
  return(list(
    chain_state = state,
    current_pdf = logpdf
  ))
}

coupledMH_kernel <-
  function(state1, state2, logpdf1, logpdf2, iteration) {
    proposal_value <-
      unbiasedmcmc::rnorm_max_coupling(state1, state2, sd_proposal, sd_proposal)
    proposal1 <- proposal_value$xy[1]
    proposal2 <- proposal_value$xy[2]
    proposal_pdf1 <- logtarget(proposal1)
    proposal_pdf2 <- logtarget(proposal2)
    logu <- log(runif(1))
    if (is.finite(proposal_pdf1)) {
      if (logu < (proposal_pdf1 - logpdf1)) {
        state1 <- proposal1
        logpdf1 <- proposal_pdf1
      }
    }
    if (is.finite(proposal_pdf2)) {
      if (logu < (proposal_pdf2 - logpdf2)) {
        state2 <- proposal2
        logpdf2 <- proposal_pdf2
      }
    }
    return(list(
      chain_state1 = state1,
      chain_state2 = state2,
      current_pdf1 = logpdf1,
      current_pdf2 = logpdf2
    ))
  }

# L values were chosen adaptively based on algorithm in paper.
repeats <- 10000
max_iterations <- 1000

pad0 <- function(x, L) {
  return(c(x, rep(0, L - length(x))))
}

# Run for a given IPM bounding function, storing results in a dataframe.
run_for_ipm_bound <- function(ipm_bound_fn) {
  bound_df_list <- list()
  for (lag in c(1, 150)) {
    bounds <- foreach(i = 1:repeats, .combine = rbind) %dorng% {
      res <-
        simulate_ipm_bounds(
          singleMH_kernel,
          coupledMH_kernel,
          sample_initial,
          ipm_bound_fn,
          max_iterations = max_iterations,
          L = lag
        )
      pad0(res$bounds, max_iterations)
    }
    ub_estimates <- colMeans(bounds)
    bound_df_list[[length(bound_df_list) + 1]] <-
      data.frame(
        lag = lag,
        t = 1:length(ub_estimates),
        bound = ub_estimates
      )
  }
  bound_df <- bind_rows(bound_df_list)
  return(bound_df)
}

tv_bound_df <- run_for_ipm_bound(tv_ipm_bound)
wasserstein_bound_df <- run_for_ipm_bound(wasserstein_ipm_bound)

########
# Plots.
########

# ggridges traceplot of approximate densities.
trace_plot_df <- data.frame(cbind(c(1:chain_length), t(all_chains)))
trace_plot_df <- melt(trace_plot_df, id.vars = "X1")
colnames(trace_plot_df) <- c('t', 'chain', 'value')
mcmc_traceplot <-
  ggplot(data = trace_plot_df, aes(x = value, y = factor(2*((t-1)%/%2)))) +
  geom_density_ridges(scale = 20, size = 0.00001, rel_min_height = 0.0005) +
  ylim(0, 160) +
  scale_x_continuous() +
  scale_y_discrete(expand = c(0,0), breaks = c(0, 50, 100, 150)) +
  xlab("x") +
  ylab("t") +
  xlim(-5, 11) +
  theme_grey(base_size = 14) +
  coord_flip()

# Total variation plot.
tv_df <- rbind(
  tv_bound_df, 
  data.frame(t = 1:chain_length - 1, lag = "exact", bound = exact_tvs)
)
g_tv <-
  ggplot(data = tv_df, aes(x = t, y = bound, linetype = lag)) +
  geom_line() +
  scale_x_continuous(expand = c(0,0), breaks = c(0, 50, 100, 150),
                     limits = c(0, 160)) +
  ylim(0, 1.5) +
  labs(y = TeX("d_{TV}")) +
  labs(x = "iteration") +
  labs(linetype = "Lag") +
  scale_linetype_manual(
    values = c(3:1),
    breaks = c("1", "150", "exact"),
    labels = c("L=1", "L=150", "Exact TV")
  ) +
  theme_grey(base_size = 14) +
  theme(legend.position = "bottom", legend.box = "horizontal")

# Wasserstein plot.
wasserstein_df <- rbind(
  wasserstein_bound_df,
  data.frame(t=1:chain_length - 1, lag = "exact", bound = exact_wasserstein)
)
g_wasserstein <-
  ggplot(data = wasserstein_df, aes(x = t, y = bound, linetype = lag)) + 
  geom_line() +
  scale_x_continuous(expand = c(0,0), breaks = c(0, 50, 100, 150),
                     limits = c(0, 160)) +
  ylim(0, 12) +
  labs(y = TeX("d_{W}")) +
  labs(x = "iteration") +
  labs(linetype = "Lag") +
  scale_linetype_manual(
    values = c(3:1),
    breaks = c("1", "150", "exact"),
    labels = c("L=1", "L=150", "Exact Wasserstein")
  ) +
  theme_grey(base_size = 14) +
  theme(legend.position = "bottom", legend.box = "horizontal")

# Combined plot legend
shared_legend <-
  cowplot::get_legend(g_wasserstein +
                        scale_linetype_manual(values = c(3:1),
                                              breaks = c("1", "150", "exact"),
                                              labels = c("L=1", "L=150",
                                                         "Exact")) +
                        theme_grey(base_size = 14) +
                        theme(legend.position = "bottom",
                              legend.box = "horizontal"))

# Vertical combined plot
tv_bound_pic_1 <- grid.arrange(mcmc_traceplot + theme(legend.position="none"),
             g_tv + theme(legend.position="none"),
             g_wasserstein + theme(legend.position="none"),
             shared_legend,
             nrow=4, ncol=1, 
             layout_matrix = cbind(c(1,2,3,4)),
             heights = c(0.31, 0.31, 0.31, 0.07))

ggsave(filename = "bound_standard_normal_vertical.pdf", plot = tv_bound_pic_1, width = 4, height = 5)


# Horizontal combined plot
tv_bound_pic_2 <- grid.arrange(mcmc_traceplot + theme(legend.position="none"),
                               g_tv + theme(legend.position="none"),
                               g_wasserstein + theme(legend.position="none"),
                               shared_legend,
                               nrow=2, ncol=3, 
                               layout_matrix = rbind(c(1,2,3), c(4,4,4)),
                               widths=c(1/3, 1/3, 1/3), heights = c(2.5, 0.2))
ggsave(filename = "bound_standard_normal_horizontal.pdf", plot = tv_bound_pic_2, width = 8, height = 3)

