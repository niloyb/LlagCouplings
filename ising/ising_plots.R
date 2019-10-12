library(latex2exp)
library(ggplot2)

source("tv_wasserstein_bounds.R")

ssg_lag <- 1e6
pt_lag <- 2e4

## TV bounds for SSG and PT at beta = 0.46

ssg_meeting_times <- read.csv("ising_ssg_meeting_times.csv", header = TRUE, sep=",")
pt_meeting_times <- read.csv("ising_pt_meeting_times.csv", header = TRUE, sep=",")

# Use only iteration numbers beginning with a single digit followed by zeros.
iterations <- c(1)
inc <- 1
while (iterations[length(iterations)] < 1e6) {
  iterations <- c(iterations, iterations[length(iterations)] + inc)
  if (length(iterations) > 1 && (((length(iterations) - 1) %% 9) == 0)) {
    inc <- inc * 10
  }
}

ssg_tv_bounds <- sapply(iterations, function(x) mean(tv_upper_bound_estimates(ssg_meeting_times$V1, ssg_lag, t = x)))
pt_tv_bounds <- sapply(iterations, function(x) mean(tv_upper_bound_estimates(pt_meeting_times$V1, pt_lag, t = x)))

ssg.df <- data.frame(x = iterations, y = ssg_tv_bounds, alg="ssg")
pt.df <- data.frame(x = iterations, y = pt_tv_bounds, alg="pt")

df <- rbind(ssg.df, pt.df)

g <- ggplot(df, aes(x = x, y = y, linetype = alg)) + 
  geom_line() +
  scale_x_log10(breaks = 10 ^ (1:6)) +
  xlab("iterations") +
  theme_grey(base_size = 14) + 
  theme(legend.position="bottom", legend.box = "horizontal") + labs(linetype="") +
  scale_linetype_manual(values = c(1, 2),
                        breaks=c("ssg", "pt"),
                        labels=c("SSG", "PT")) +
  labs(y = TeX("d_{TV}"))
g
