L-Lag Couplings
===============

These scripts reproduce the results of article "Estimating Convergence of Markov chains with L-Lag Couplings", by Niloy Biswas, Pierre E. Jacob and Paul Vanetti. In Advances in Neural Information Processing Systems 32 (2019), 7389--7399 (https://papers.nips.cc/paper/8958-estimating-convergence-of-markov-chains-with-l-lag-couplings). 


## Tutorial
An R Markdown tutorial is available: https://niloyb.github.io/LlagCouplings/

## References
The unbiasedmcmc package is based on "Unbiased Markov chain Monte Carlo with couplings", by Pierre E. Jacob, John O'Leary, Yves F. Atchade. https://arxiv.org/abs/1708.03625

The debiasedhmc package is based on "Unbiased Hamiltonian Monte Carlo with couplings", by Jeremy Heng and Pierre E. Jacob. https://academic.oup.com/biomet/article/106/2/287/5366709

The Polya-Gamma sampler is taken from the package BayesLogit of Nick Polson, James Scott, and Jesse Windle. https://cran.r-project.org/package=BayesLogit

## Erratum
There was a bug in the code used to produce Figure 5 of the above article. The choice of step-size for the MALA and ULA chains was not set in the code as it was described in the paper. The corresponding script has now been corrected. The folder *ula_versus_mala* now contains the corrected script, updated simulation results and a note with details.