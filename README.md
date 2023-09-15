# GARP

R codes for graph-aligned random partition (GARP) model.

**Authors**: [Giovanni Rebaudo](https://giovannirebaudo.github.io) and [Peter Müller](https://web.ma.utexas.edu/users/pmueller).

#### Overview 
This repository is associated with the article [Rebaudo, G. and Müller, P. (2023). Graph-aligned Random Partition Model (GARP). Submitted.](https://arxiv.org/abs/2306.08485)
The key contribution of the paper is outlined below.
 
> [...] Motivated by single-cell RNA applications we develop a novel dependent mixture model to jointly perform cluster analysis and align the clusters on a graph.
Our flexible graph-aligned random partition model (GARP) exploits Gibbs-type priors as building blocks, allowing us to derive analytical results on the graph-aligned random partition's probability mass function (pmf).
We derive a generalization of the Chinese restaurant process from the pmf and a related efficient and neat MCMC algorithm to perform Bayesian inference. 

This repository provides codes to replicate the results in Rebaudo and Müller (2023). Graph-aligned Random Partition Model (GARP).

In particular, we provide the `R` code to implement the MCMC to perform posterior inference under the GARP model.

The repository contains the following:

1. `GARP_main.R` code to reproduce the main results in the article;
2. `GARP_fcts.R` functions needed to run the main code;
3. `Data-and-Results` folder with data and results of the analyses.

#### Questions or bugs
For bug reporting purposes, e-mail [Giovanni Rebaudo](https://giovannirebaudo.github.io) (giovanni.rebaudo@gmail.com)

#### Citation
Please cite the following publication if you use this repository in your research: [Rebaudo, G. and Müller, P. (2023). Graph-aligned Random Partition Model (GARP). Submitted.](https://arxiv.org/abs/2306.08485)




