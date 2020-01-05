### Unadjusted and Metropolis Adjusted Langevin Algorithm
===============

These scripts reproduce the ULA (Unadjusted Langevin Algorithm) versus MALA (Metropolis Adjusted Langevin Algorithm) example. That is, obtain meeting times for the ULA and MALA samplers targetting a d-dimensional multivariate gaussian with mean zero and covariance matrix [Sigma]{i,j}=0.5^|i-j| for 1<= i, j <=d, for d = 50, 100, 200, 300, 400, 500, 600, 800, 1000. The standard deviation of the ULA and MALA proposals are 0.1d^{-1/6} and d^{-1/6} respectively.

* `ula_versus_mala.R` contains all the functions and generates meeting times for ULA and MALA for varying dimensions. It requires the packages `doParallel`, `doRNG` and `unbiasedmcmc`, and the source file `coupling.R`.
* `ula_versus_mala_plots.R` plots the figure from the paper. It requires the packages `ggplot2` and `latex2exp`, and the source file `tv_wasserstein_bounds.R`.

