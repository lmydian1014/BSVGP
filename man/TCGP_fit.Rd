\name{TCGP_fit}
\alias{TCGP_fit}
\title{
Bayesian fitting of the Thresholded Correlation Gaussian Process.
}

\usage{
TCGP_fit(
  grids, 
  T, 
  V, 
  n, 
  L, 
  Xmat, 
  lambda, 
  tau_1_sq_init,
  tau_2_sq_init, 
  c_init, 
  e_pos_init, 
  e_neg_init, 
  Y_pos, 
  Y_neg, 
  alpha1, 
  beta1, 
  alpha2, 
  beta2, 
  thres_init, 
  rinvgamma)
}

\arguments{
  \item{girds}{A matrix of real numbers as grid points where rows are observations and columns are coordinates.}
  \item{T}{A integer number to specify the number of iteractions in MCMC sampling.}
  \item{V}{A integer number to specify the number of voxels.}
  \item{n}{A integer number to specify the number of subjects.}
  \item{L}{A integer number to specify the number of basis functions.}
  \item{Xmat}{A matrix represents the basis functions evaluated at the grid points, where rows are observations and columns are the basis functions.}
  \item{lambda}{A vector represents the eigen values.}
  \item{tau_1_sq_init}{A vector of length V specifies the initial value of tau_1_sq.}
  \item{tau_2_sq_init}{A vector of length V specifies the initial value of tau_2_sq.}
  \item{c_init}{A vector of length L specifies the initial value of c.}
  \item{e_pos_init}{A matrix of dimension L*n specifies the initial value of e_pos.}
  \item{e_neg_init}{A matrix of dimension L*n specifies the initial value of e_neg.}
  \item{Y_pos}{A matrix of dimension V*n specifies the transformed image, where Y_pos = (Y_1 + Y_2)/2.}
  \item{Y_neg}{A matrix of dimension V*n specifies the transformed image, where Y_neg = (Y_1 - Y_2)/2.}
  \item{alpha1}{A positive real number specify the shape parameter in the inverse Gamma distribution, which is the prior of tau_1_sq.}
  \item{beta1}{A positive real number specify the rate parameter in the inverse Gamma distribution, which is the prior of tau_1_sq.}
  \item{alpha2}{A positive real number specify the shape parameter in the inverse Gamma distribution, which is the prior of tau_2_sq.}
  \item{beta2}{A positive real number specify the rate parameter in the inverse Gamma distribution, which is the prior of tau_2_sq.}
  \item{thres_init}{A non-negative real number specify the initial value of the thresholding parameter thres.}
  \item{rinvgamma}{R function rinvgamma.}

}


\value{
  \item{gibbs_c }{A matrix of dimension L*T represents the posterior samples of c for each iteration.}
  \item{gibbs_tau_1_sq }{A matrix of dimension V*T represents the posterior samples of tau_1_sq for each iteration.}
  \item{gibbs_tau_2_sq }{A matrix of dimension V*T represents the posterior samples of tau_2_sq for each iteration.}
  \item{temp_logL }{A vector of length T record the loglikelihood for each iteration.}
  \item{gibbs_thres }{A vector of length T represents the posterior samples of the thresholding parameter for each iteration.}
}


\author{
Moyan Li <moyanli@umich.edu>
}

\examples{
  dat = gen_data_design(n=50, d = 2, num_grids = 32, grids_lim = c(0,1), poly_degree = 32, a = 0.1, b = 1, pos_radius = 0.1, neg_radius = 0.1, pos_mag = 0.75, neg_mag = 0.85)
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
  T = 800
  chain = TCGP_fit(grids, T, V, n, L, Xmat, lambda, tau_1_sq_init,tau_2_sq_init, c_init, e_pos_init, e_neg_init, Y_pos, Y_neg, 1, 1, 1, 1, 0, rinvgamma)
}
