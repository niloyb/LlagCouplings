### Ising model 
rm(list = ls())
set.seed(21)
#
library(doRNG)
library(doParallel)
library(latex2exp)
library(ggplot2)
registerDoParallel(cores = detectCores()-2)
library(Rcpp)

## 
source("tv_wasserstein_bounds.R")

## TV bounds for SSG and PT at beta = 0.46
load(file = "ising.tvbounds.pt.RData")
load(file = "ising.tvbounds.ssg.RData")

### now comparison of the two samplers (SSG in full lines, PT in dashed lines)
## the following takes a minute to run, due to the number of bounds (1e6)
## this can probably be sped up without noticeable difference in the plot
## e.g. by subsampling the iterations
## in the meantime we can save the results for faster re-plotting

iterations <- 1:1e6 
# ssg_tv_bounds <- sapply(iterations, function(x) mean(tv_upper_bound_estimates(ssg_meetingtimes2, ssg_lag, t = x)))
# pt_tv_bounds <- sapply(iterations, function(x) mean(tv_upper_bound_estimates(pt_meetingtimes2, pt_lag, t = x)))
# save(ssg_tv_bounds, pt_tv_bounds, file = "ising.bothbounds.RData")
load(file = "ising.bothbounds.RData")

plot(ssg_tv_bounds, type = "l", xlab = "iteration", ylab = "tv bounds", ylim = c(0,1.1), log = "x")
lines(x = iterations, y = pt_tv_bounds, lty = 2)
abline(h = 0)

ssg.df <- data.frame(x = iterations, y = ssg_tv_bounds)
pt.df <- data.frame(x = iterations, y = pt_tv_bounds)

g <- ggplot(ssg.df, aes(x = x, y = y)) + geom_line(aes(linetype = "SSG")) + scale_x_log10(breaks = c(1e1, 1e2, 1e3, 1e4, 1e5, 1e6)) +
  geom_line(data=pt.df, aes(linetype = "PT")) + xlab("iterations")
g <- g + theme_grey(base_size = 14) + 
  theme(legend.position="bottom", legend.box = "horizontal") + labs(linetype="") + 
  scale_linetype_manual(values = c(1,2),
                        breaks=c("SSG", "PT"),
                        labels=c("SSG", "PT"))
g <- g + labs(y = TeX("d_{TV}"))
g
ggsave(filename = "/Users/niloybiswas/Dropbox/total_variation_coupling/NeurIPS2019 write up/images/ising_ssg_versus_pt_2.pdf", plot = g, width = 8, height = 4)


# 
# ## Mixing time of single site Gibbs (SSG) as a function of beta
# ## the mixing time is defined as the first iteration for which 
# ## the estimate TV bound is less than epsilon = 0.1 
# load(file = "ising.mixingtime.vs.beta.RData")
# # based on NREP replications for each beta in a grid
# plot(betas_, mixingtimes, type = "b", xlab = "beta", ylab = "mixing time")
# 
# ## and same thing for parallel tempering
# load(file = "ising.pt.mixingtime.vs.beta.RData")
# plot(betas_, mixingtimes, type = "l", xlab = "beta", ylab = "mixing time", log = "y", xlim = c(0.3, 0.5))
# lines(ptbetas, ptmixingtimes, lty = 2)
