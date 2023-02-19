# Bayesian Inferences on Spatially-Varying Correlations via Gaussian Processes
An R package for performing the Bayesian Inferences on Spatially-Varying Correlations via Gaussian Processes

- Install and load the package
  ```
  # Please make sure you have installed the following packages
  # truncnorm, MASS, invgamma, BayesGPfit, lattice, sn, coda, mcmcplots, mcmc, Rcpp, plgp, matrixStats, igraph
  devtools::install_github("lmydian1014/BSV-GP")
  ```
- Simulate the image
 ```
  dat = gen_data(n = 50, d = 2, num_grids = 8, grids_lim = c(-1,1), random = FALSE, 
    poly_degree = 32, a = 0.01, b = 10, center = NULL, rate = NULL, max_range = 6, thres = 1)
  
  dat = gen_data_design(n = 50, d = 2, num_grids = 12, grids_lim = c(0,1), poly_degree = 32, a = 0.01, b = 10, 
    pos_radius = pos_radius, neg_radius = neg_radius,
    pos_mag = 0.75, neg_mag = 0.85)
 ```
- Bayesian Inferences on Spatially-Varying Correlations
