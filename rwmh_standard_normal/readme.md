### Standard normal target
===============

These scripts reproduce the Univariate Normal RWMH example. That is, obtain meeting times for the RWMH sampler with step-size 0.5 targeting N(0,1) and plot the corresponding TV and 1-Wasserstein bounds.

* `standard_normal_RWMH.R` both runs the simulations and creates the figures. It requires the packages `doParallel`, `doRNG`, `dplyr`, `ggplot2`, `ggridges`, `gridExtra`, `latex2exp` and `reshape2`, and the source files `coupling.R` and `tv_wasserstein_bounds.R`.

