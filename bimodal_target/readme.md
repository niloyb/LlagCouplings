### Bimodal target
===============

These scripts reproduce the Bimodal target example. That is, a histogram of the 500th marginal of a Random-Walk Metropolis--Hastings Markov chain with: initial distribution N(10,1) and step-size 1 targeting 0.5N(-4,1)+0.5N(+4,1), the meeting times of the coupled chain and corresponding TV bounds.

* `bimodal_target.R` contains all the functions and generates the data. It requires the packages `doParallel` and `doRNG`.
* `bimodal_target_plots.R` plots the figure in the paper. It requires packages `dplyr`, `ggplot2`, `latex2exp`, `purrr` and `cowplot`, and the source file `tv_wasserstein_bounds.R`.

