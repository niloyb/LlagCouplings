### Logistic regression
===============

These scripts reproduce the reproduce the Logistic Regression example. That is, obtain meeting times for (i) a coupled Polya-Gamma sampler, (ii) coupled HMC samplers with leapfrog integration step size 0.025, leapfrog number of steps 4,5,6,7, and RWMH kernel step size 0.001, and probability of selecting RWMH kernel 0.05. Both Polya-Gamma and HMC samplers target a Bayesian logistic gression posterior. 

* `hmc_logistic_regression.R` contains all the functions and generates meeting times for HMC. It requires the packages `doParallel`,
`unbiasedmcmc` and `debiasedhmc`, and the source file `coupling.R`.
* `polya_gamma.R` contains all the functions and generates meeting times for the Polya-Gamma sampler. It requires the packages `doParallel`, `unbiasedmcmc` and `BayesLogit`, and the source file `coupling.R`.
* `logistic_regression_plots.R` plots the figures from the paper. It requires the packages `latex2exp` and `tidyverse`, and the source file `tv_wasserstein_bounds.R`.

