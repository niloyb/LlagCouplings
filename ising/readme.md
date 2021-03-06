### Ising model
===============

These scripts reproduce the Ising model example. That is, obtain meeting times for the single site Gibbs (SSG) with a lag of 10^6 at beta = 0.46 and meeting times for parallel tempering (PT) with a lag of 2x10^4 with 12 chains on betas = 0.3,...,0.46 (equispaced).

* `ising.cpp` contains functions coded using Rcpp to perform Gibbs steps and couplings thereof on an Ising model.
* `ising_functions.R` contains R functions, some of which call the Rcpp functions defined in `ising.cpp`.
* `ising_run.R` runs the simulations. It produces two files: ising.tvbounds.ssg.RData and ising.tvbounds.pt.RData. The time to complete depends on the number of cores; on a MacBook pro from 2015 and using 7 cores, it took 6 hours (for 500 independent replicates). It requires the packages `doParallel` and `doRNG`, and the source files `coupling.R` and `ising_functions.R`.
* `ising_plots.R` loads the above files, which contain the meeting times, and creates the figure of the paper. Note that since there are 10^6 bounds to compute this takes a minute or so. It requires the packages `ggplot2`, `latex2exp`, and the source file `tv_wasserstein_bounds.R`.