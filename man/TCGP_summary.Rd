\name{TCGP_summary}
\alias{TCGP_summary}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 The summary of the posterior sampling results of the Thresholded Correlation Gaussian Process.
}

\usage{
TCGP_summary(T, dat, chain, burn_in, grids, Xmat, thres)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{T}{A integer number to specify the number of iteractions in MCMC sampling.}
  \item{dat}{A list of data including grids, Y_1 and Y_2.}
  \item{chain}{A list of the posterior sampling results obtained from TCGP_fit().}
  \item{burn_in}{An integer number to specify the burn-in number.}
  \item{grids}{A matrix of real numbers as grid points where rows are observations and columns are coordinates.}
  \item{Xmat}{A matrix real numbers for the basis functions evaluated at the grid points, where rows are observations and columns are the basis functions.}
}


\value{
    \item{cor_type}{A vector with length V specifies the correlation type of each voxel. cor_type = 1 and -1 represents the significant positive or negative correlation while cor_type = 0 represents there is no significant correlation.}
    \item{rho_hat}{A vector with length V specifies the mean correlation of each voxel, which is the mean of the posterior samples of the correlation.}
    \item{sensi_pos}{The sensitivity of the positive correlated region.}
    \item{sensi_neg}{The sensitivity of the negative correlated region.}
    \item{speci_pos}{The specificity of the positive correlated region.}
    \item{speci_neg}{The specificity of the negative correlated region.}
    \item{FDR_pos}{The False Discovery Rate of the positive correlated region.}
    \item{FDR_neg}{The False Discovery Rate of the negative correlated region.}
    \item{prob_pos}{The probability map of the positive correlated region.}
    \item{prob_neg}{The probability map of the negative correlated region.}
    \item{prob_0}{The probability map of the non-correlated region.}
}

\author{
Moyan Li <moyanli@umich.edu>
}

\examples{
### generate data
  dat = gen_data_design(n=50, d = 2, num_grids = 32, grids_lim = c(0,1), poly_degree = 32, a = 0.1, b = 1, pos_radius = 0.1, neg_radius = 0.1,   pos_mag = 0.75, neg_mag = 0.85)
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
  
### obtain the MCMC chain
  T = 800
  chain = TCGP_fit(grids, T, V, n, L, Xmat, lambda, tau_1_sq_init,tau_2_sq_init, c_init, e_pos_init, e_neg_init, Y_pos, Y_neg, 1, 1, 1, 1, 0, rinvgamma)
  
### summary results
  res = TCGP_summary(T = T, chain = chain, dat = dat, burn_in = 0.2*T, grids, Xmat)
  
### plot 
  grid.panel = function(...) {
    panel.levelplot(...)
    panel.abline(h = seq(0.25, 0.75, length = 3), v = seq(0.25, 0.75, length = 3), lty = 3, col = gray(0.5))
  }
  fig = fourfigs.levelplot(
    res$prob_pos, res$prob_neg, dat$rho, res$cor_type,
    grids[, 1],grids[, 2],
    titles = c("Pos.Cor.Prob", "Neg.Cor.Prob", "True correlation", "Selection"),
    layout = c(2, 2), panel = grid.panel)
  }
