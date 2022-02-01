### Unadjusted and Metropolis Adjusted Langevin Algorithm
===============

These scripts reproduce the ULA (Unadjusted Langevin Algorithm) versus MALA (Metropolis Adjusted Langevin Algorithm) example. 

* `ula_versus_mala.R` contains all the functions and generates meeting times for ULA and MALA for varying dimensions, and plots the figures. It requires the packages `doParallel`, `doRNG`, `ggplot2` and `latex2exp`, and the source files `coupling.R`, `mvnorm.R`, `mvnorm_couplings.R`, `mvnorm.cpp` and `tv_wasserstein_bounds.R`.

* `fig5_erratum.pdf' contains details about a bug on the previous version of the code, which has now been corrected.