# Bayesian Inferences on Spatially-Varying Correlations via the Thresholded Correlation Gaussian Processes
An R package for performing the Bayesian Inferences on Spatially-Varying Correlations via the Thresholded Correlation Gaussian Processes. 

### Install and load the package
```
list.of.packages = c("truncnorm", "MASS", "invgamma", "BayesGPfit", "lattice", "sn", "coda", "mcmcplots", "mcmc","Rcpp", "plgp", "matrixStats", "igraph")
new.packages = list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)>0) install.packages(new.packages)
  
devtools::install_github("lmydian1014/BSVGP")
library(BSVGP)

```
### Simulate the image
 ```
dat = gen_data_design(n=50, d = 2, num_grids = 32, grids_lim = c(0,1), poly_degree = 32, a = 0.1, b = 1, pos_radius = 0.1, neg_radius = 0.1, pos_mag = 0.75, neg_mag = 0.85)
 ```
### Bayesian Inferences on Spatially-Varying Correlations

#### model representation and generate initial values for parameters of interest
```
n = 50
Y_1 = dat$Y_1
Y_2 = dat$Y_2 
Y_pos = (Y_1 + Y_2)/2
Y_neg = (Y_1 - Y_2)/2
grids = dat$x
eig = approx_eigen(n = n, grids, 64, l = 0.05)
Xmat = eig$Xmat
lambda = eig$lambda
L = length(lambda)
V = nrow(grids)
e_pos_init = matrix(NA, nrow = length(lambda), ncol = n) 
e_neg_init = matrix(NA, nrow = length(lambda), ncol = n) 
for(i in c(1:n)){
    e_pos_init[,i] = rnorm(length(lambda), 0, sqrt(lambda))
    e_neg_init[,i] = rnorm(length(lambda), 0, sqrt(lambda))
}
c_init = rnorm(length(lambda), 0, sqrt(lambda))
tau_1_sq_init = rinvgamma(V, 3, 0.1)
tau_2_sq_init = rinvgamma(V, 3, 0.1)
```
#### Run BSV-GP model
```
T = 800
chain = sample_gibbs_cpp(grids, T, V, n, L, Xmat, lambda, tau_1_sq_init,tau_2_sq_init, c_init, e_pos_init, e_neg_init, Y_pos, Y_neg, 1, 1, 1, 1, 0, rinvgamma)
```
#### Analysis the MCMC chain and perform selection
```
res_gibbs = analysis_chain(T = T, chain = chain, dat = dat, burn_in = 0.2*T, grids, Xmat, thres)
```
#### Plot the selection results and the inclusion probability map
```
grid.panel = function(...) {
    panel.levelplot(...)
    panel.abline(h = seq(0.25, 0.75, length = 3), v = seq(0.25, 0.75, length = 3), lty = 3, col = gray(0.5))
}

fig = fourfigs.levelplot(
    res_gibbs$prob_pos, res_gibbs$prob_neg, dat$rho, res_gibbs$cor_type,
    grids[, 1],grids[, 2],
    titles = c("Pos.Cor.Prob", "Neg.Cor.Prob", "True correlation", "Selection"),
    layout = c(2, 2), panel = grid.panel)
```

![alt text](https://github.com/lmydian1014/BSVGP/blob/main/fig.png)

