#'@title The summary of the posterior sampling results of the Thresholded Correlation Gaussian Process.

#'@param T A integer number to specify the number of iteractions in MCMC sampling.
#'@param dat A list of data including grids, Y_1 and Y_2.
#'@param chain A list of the posterior sampling results obtained from TCGP_fit().
#'@param burn_in An integer number to specify the burn-in number.
#'@param grids A matrix of real numbers as grid points where rows are observations and columns are coordinates.
#'@param Xmat A matrix real numbers for the basis functions evaluated at the grid points, where rows are observations and columns are the basis functions.

#'
#'@return A list of variables including the model fitting results
#'\describe{
#'\item{cor_type}{A vector with length V specifies the correlation type of each voxel. cor_type = 1 and -1 represents the significant positive or negative correlation while cor_type = 0 represents there is no significant correlation.}
#'\item{rho_hat}{A vector with length V specifies the mean correlation of each voxel, which is the mean of the posterior samples of the correlation.}
#'\item{sensi_pos}{The sensitivity of the positive correlated region.}
#'\item{sensi_neg}{The sensitivity of the negative correlated region.}
#'\item{speci_pos}{The specificity of the positive correlated region.}
#'\item{speci_neg}{The specificity of the negative correlated region.}
#'\item{FDR_pos}{The False Discovery Rate of the positive correlated region.}
#'\item{FDR_neg}{The False Discovery Rate of the negative correlated region.}
#'\item{prob_pos}{The probability map of the positive correlated region.}
#'\item{prob_neg}{The probability map of the negative correlated region.}
#'\item{prob_0}{The probability map of the non-correlated region.}
#'}
#'
#'@author Moyan Li <moyanli@umich.edu>
#'
#'@examples
#'\examples{
#'  dat = gen_data_design(n=50, d = 2, num_grids = 32, grids_lim = c(0,1), poly_degree = 32, a = 0.1, b = 1, pos_radius = 0.1, neg_radius = 0.1, pos_mag = 0.75, neg_mag = 0.85)
#'  n = 50
#'  Y_1 = dat$Y_1
#'  Y_2 = dat$Y_2 
#'  Y_pos = (Y_1 + Y_2)/2
#'  Y_neg = (Y_1 - Y_2)/2
#'  grids = dat$x
#'  eig = approx_eigen(n = n, grids, 64, l = 0.05)
#'  Xmat = eig$Xmat
#'  lambda = eig$lambda
#'  L = length(lambda)
#'  V = nrow(grids)
#'  e_pos_init = matrix(NA, nrow = length(lambda), ncol = n) 
#'  e_neg_init = matrix(NA, nrow = length(lambda), ncol = n) 
#'  for(i in c(1:n)){
#'      e_pos_init[,i] = rnorm(length(lambda), 0, sqrt(lambda))
#'      e_neg_init[,i] = rnorm(length(lambda), 0, sqrt(lambda))
#'  }
#'  c_init = rnorm(length(lambda), 0, sqrt(lambda))
#'  tau_1_sq_init = rinvgamma(V, 3, 0.1)
#'  tau_2_sq_init = rinvgamma(V, 3, 0.1)
#'  T = 800
#'  chain = TCGP_fit(grids, T, V, n, L, Xmat, lambda, tau_1_sq_init,tau_2_sq_init, c_init, e_pos_init, e_neg_init, Y_pos, Y_neg, 1, 1, 1, 1, 0, rinvgamma)
#'}
#'  res = TCGP_summary(T = T, chain = chain, dat = dat, burn_in = 0.2*T, grids, Xmat)
#'
#'@export

TCGP_summary = function(T, dat, chain, burn_in, grids, Xmat){
	gibbs_c = chain$gibbs_c
	gibbs_tau_1_sq = chain$gibbs_tau_1_sq
	gibbs_tau_xi_sq = 1
	gibbs_tau_2_sq = chain$gibbs_tau_2_sq
	gibbs_thres = chain$gibbs_thres
	temp_logL = chain$temp_logL
	rho_hat = matrix(nrow = nrow(grids), ncol = T)
	rho_mean = rep(0, nrow(grids))
	xi_hat_m = matrix(nrow = nrow(grids), ncol = T)
	sigma_pos_sq_hat_m = matrix(nrow = nrow(grids), ncol = T)
	sigma_neg_sq_hat_m = matrix(nrow = nrow(grids), ncol = T)
	for(i in c(1:T)){
		xi_hat_m[,i] = Xmat %*% gibbs_c[,i]
		sigma_pos_sq_hat_m[,i] = ifelse(xi_hat_m[,i] > gibbs_thres[i],xi_hat_m[,i], 0)
        sigma_neg_sq_hat_m[,i] = ifelse(xi_hat_m[,i] < -gibbs_thres[i],-xi_hat_m[,i], 0)
		numer = sigma_pos_sq_hat_m[,i] - sigma_neg_sq_hat_m[,i]
		denor1 = sqrt(sigma_pos_sq_hat_m[,i] + sigma_neg_sq_hat_m[,i] + gibbs_tau_1_sq[,i])
		denor2 = sqrt(sigma_pos_sq_hat_m[,i] + sigma_neg_sq_hat_m[,i] + gibbs_tau_2_sq[,i])
		rho_hat[,i]=  numer/(denor1 * denor2)
	}
	rho_hat = rho_hat[,burn_in:T]
	rho_true = dat$rho 
	rho_true[dat$rho > 0] = 1
	rho_true[dat$rho < 0] = -1
	prob_pos = c(); prob_neg = c(); prob_0 = c()
	sensi_pos_num = sensi_neg_num = speci_pos_num = speci_neg_num = 0
	FDR_pos_num = FDR_neg_num = 0
	for(v in c(1:nrow(grids))){
		num_0 = sum(rho_hat[v,] == 0)
		num_pos = sum(rho_hat[v,] > 0)
		num_neg = sum(rho_hat[v,] < 0)
		prob_pos = append(prob_pos, num_pos/(T - burn_in + 1))
		prob_neg = append(prob_neg, num_neg/(T - burn_in + 1))
		prob_0 = append(prob_0, num_0/(T - burn_in + 1))
	  
		if(num_0 > (T - burn_in + 1)/3){
			next;
		}
	    else if(num_pos > num_neg){
	    	rho_mean[v] = 1
	    }
	    else{
	    	rho_mean[v] = -1
	    }
	}
	for(v in c(1:nrow(grids))){
		if(rho_true[v] > 0){
			if(rho_mean[v] > 0){
				sensi_pos_num = sensi_pos_num + 1
				speci_neg_num = speci_neg_num + 1
			}
			if(rho_mean[v] == 0){
				speci_neg_num = speci_neg_num + 1
			}
			# if(rho_mean[v] < 0){
			# 	FDR_neg_num = FDR_neg_num + 1
			# }
		}
		if(rho_true[v] == 0){
			if(rho_mean[v] == 0){
				speci_pos_num = speci_pos_num + 1
				speci_neg_num = speci_neg_num + 1
			}
			if(rho_mean[v] > 0){
				speci_neg_num = speci_neg_num + 1
				FDR_pos_num = FDR_pos_num + 1
			}
			if(rho_mean[v] < 0){
				speci_pos_num = speci_pos_num + 1
				FDR_neg_num = FDR_neg_num + 1
			}
			
		}
		if(rho_true[v] < 0){
			if(rho_mean[v] < 0){
				sensi_neg_num = sensi_neg_num + 1
				speci_pos_num = speci_pos_num + 1
			}
			if(rho_mean[v] == 0){
				speci_pos_num = speci_pos_num + 1
			}
			# if(rho_mean[v] > 0){
			# 	FDR_pos_num = FDR_pos_num + 1
			# }
		}
	}
	sensi_pos = sensi_pos_num/sum(rho_true == 1)
	sensi_neg = sensi_neg_num/sum(rho_true == -1)
	speci_pos = speci_pos_num/(nrow(grids) - sum(rho_true == 1))
	speci_neg = speci_neg_num/(nrow(grids) - sum(rho_true == -1))
	FDR_pos = FDR_pos_num /(sensi_pos_num + FDR_pos_num)
	FDR_neg = FDR_neg_num /(sensi_neg_num + FDR_neg_num)

	tau_1_sq_hat = rowSums(gibbs_tau_1_sq[,burn_in:T])/(T-burn_in +1)
	tau_2_sq_hat = rowSums(gibbs_tau_2_sq[,burn_in:T])/(T-burn_in +1)

	rho2 = rowSums(rho_hat)/(T-burn_in+1)
	xi_hat = rowSums(xi_hat_m[, c(burn_in:T)])/(T-burn_in +1)
	sigma_pos_sq_hat = rowSums(sigma_pos_sq_hat_m[, c(burn_in:T)])/(T-burn_in +1)
	sigma_neg_sq_hat = rowSums(sigma_neg_sq_hat_m[, c(burn_in:T)])/(T-burn_in +1)

	thres_xi = sigma_pos_sq_hat - sigma_neg_sq_hat
	abs_thres_xi = abs(thres_xi)
	var_Y_1_hat = tau_1_sq_hat + abs_thres_xi
	var_Y_2_hat = tau_2_sq_hat + abs_thres_xi
	cov_hat = sigma_pos_sq_hat - sigma_neg_sq_hat
	true_cov = dat$sigma_pos_sq - dat$sigma_neg_sq

	true_thres_xi = dat$sigma_pos_sq - dat$sigma_neg_sq
	true_abs_thres_xi = abs(true_thres_xi)
	true_var_Y_1 = dat$tau_1_sq + true_abs_thres_xi
	true_var_Y_2 = dat$tau_2_sq + true_abs_thres_xi
	
	#return(list('xi_hat_m' = xi_hat_m, 'cor_type' = rho_mean, 'rho_hat' = rho_hat, 'var_Y_1_hat' = var_Y_1_hat, 'var_Y_2_hat' = var_Y_2_hat, 'cov_hat' = cov_hat, 
	#	'fig_rho' = fig_rho, 'fig_rho2' = fig_rho2,'fig_sigma_pos_sq' = fig_sigma_pos_sq, 'fig_sigma_neg_sq' = fig_sigma_neg_sq, 
	#	'fig_cov' = fig_cov, 'sensi_pos' = sensi_pos, 'sensi_neg' = sensi_neg, 'speci_pos' = speci_pos, 'speci_neg' = speci_neg, 'FDR_pos' = FDR_pos, 'FDR_neg' = FDR_neg, 'fig_tau1' = fig_tau1,
	#	'fig_tau2' = fig_tau2, 'fig_prob_map' = fig_prob_map, "prob_pos" = prob_pos, "prob_neg" = prob_neg, "prob_0" = prob_0))
	return(list('cor_type' = rho_mean, 'rho_hat' = rho_hat, 'sensi_pos' = sensi_pos, 'sensi_neg' = sensi_neg, 'speci_pos' = speci_pos, 'speci_neg' = speci_neg, 
		'FDR_pos' = FDR_pos, 'FDR_neg' = FDR_neg, "prob_pos" = prob_pos, "prob_neg" = prob_neg, "prob_0" = prob_0))

}

analysis_chain_ABCD = function(T, chain, burn_in, V, Xmat){
	gibbs_c = chain$gibbs_c
	gibbs_tau_1_sq = chain$gibbs_tau_1_sq
	gibbs_tau_xi_sq = 1
	gibbs_tau_2_sq = chain$gibbs_tau_2_sq
	gibbs_thres = chain$gibbs_thres
	#gibbs_thres = rep(thres, T)
	temp_logL = chain$temp_logL
	rho_hat = matrix(nrow = V, ncol = T)
	rho_mean = rep(0, V)
	xi_hat_m = matrix(nrow = V, ncol = T)
	sigma_pos_sq_hat_m = matrix(nrow = V, ncol = T)
	sigma_neg_sq_hat_m = matrix(nrow = V, ncol = T)
	for(i in c(1:T)){
		xi_hat_m[,i] = Xmat %*% gibbs_c[,i]
		sigma_pos_sq_hat_m[,i] = ifelse(xi_hat_m[,i] > gibbs_thres[i],xi_hat_m[,i], 0)
  		sigma_neg_sq_hat_m[,i] = ifelse(xi_hat_m[,i] < -gibbs_thres[i],-xi_hat_m[,i], 0)

		numer = sigma_pos_sq_hat_m[,i] - sigma_neg_sq_hat_m[,i]
		denor1 = sqrt(sigma_pos_sq_hat_m[,i] + sigma_neg_sq_hat_m[,i] + gibbs_tau_1_sq[,i])
		denor2 = sqrt(sigma_pos_sq_hat_m[,i] + sigma_neg_sq_hat_m[,i] + gibbs_tau_2_sq[,i])
		rho_hat[,i]=  numer/(denor1 * denor2)
	}
	rho_hat = rho_hat[,burn_in:T]
	rho_esti = rowSums(rho_hat)/(T-burn_in+1)
	prob_pos = c(); prob_neg = c(); prob_0 = c()
	for(v in c(1:V)){
		num_0 = sum(rho_hat[v,] == 0)
		num_pos = sum(rho_hat[v,] > 0)
		num_neg = sum(rho_hat[v,] < 0)
		prob_pos = append(prob_pos, num_pos/(T - burn_in + 1))
		prob_neg = append(prob_neg, num_neg/(T - burn_in + 1))
		prob_0 = append(prob_0, num_0/(T - burn_in + 1))
		if(num_0 > (T - burn_in + 1)/4){
			next;
		}
	    if(num_pos > num_neg){
	    	rho_mean[v] = 1
	    }
	    else if(num_pos < num_neg){
	    	rho_mean[v] = -1
	    }
	    else{
	    	next;
	    }
	}
	return(cbind(prob_pos, prob_neg, prob_0, rho_mean, rho_esti))
}











analysis_chain_soft = function(T, dat, chain, burn_in, grids, Xmat, thres){
	####### input a list 

	gibbs_c = chain$gibbs_c
	gibbs_tau_1_sq = chain$gibbs_tau_1_sq
	gibbs_tau_xi_sq = 1
	gibbs_tau_2_sq = chain$gibbs_tau_2_sq
	gibbs_thres = chain$gibbs_thres
	#gibbs_thres = rep(thres, T)
	temp_logL = chain$temp_logL
	rho_hat = matrix(nrow = nrow(grids), ncol = T)
	rho_mean = rep(0, nrow(grids))
	xi_hat_m = matrix(nrow = nrow(grids), ncol = T)
	sigma_pos_sq_hat_m = matrix(nrow = nrow(grids), ncol = T)
	sigma_neg_sq_hat_m = matrix(nrow = nrow(grids), ncol = T)
	xi_hat_m = Xmat %*% gibbs_c
	xi_hat = rowSums(xi_hat_m[, c(burn_in:T)])/(T-burn_in +1)
	thres = quantile(abs(xi_hat), 0.85)
	for(i in c(1:T)){

		sigma_pos_sq_hat_m[,i] = ifelse(xi_hat_m[,i] > thres, xi_hat_m[,i] - thres, 0)
  		sigma_neg_sq_hat_m[,i] = ifelse(xi_hat_m[,i] < -thres, -xi_hat_m[,i] - thres, 0)

		numer = sigma_pos_sq_hat_m[,i] - sigma_neg_sq_hat_m[,i]
		denor1 = sqrt(sigma_pos_sq_hat_m[,i] + sigma_neg_sq_hat_m[,i] + gibbs_tau_1_sq[,i])
		denor2 = sqrt(sigma_pos_sq_hat_m[,i] + sigma_neg_sq_hat_m[,i] + gibbs_tau_2_sq[,i])
		rho_hat[,i]=  numer/(denor1 * denor2)
	}
	rho_hat = rho_hat[,burn_in:T]
	rho_true = dat$rho 
	rho_true[dat$rho > 0] = 1
	rho_true[dat$rho < 0] = -1
	prob_pos = c(); prob_neg = c(); prob_0 = c()
	sensi_pos_num = sensi_neg_num = speci_pos_num = speci_neg_num = 0
	FDR_pos_num = FDR_neg_num = 0
	for(v in c(1:nrow(grids))){
		num_0 = sum(rho_hat[v,] == 0)
		num_pos = sum(rho_hat[v,] > 0)
		num_neg = sum(rho_hat[v,] < 0)
		prob_pos = append(prob_pos, num_pos/(T - burn_in + 1))
		prob_neg = append(prob_neg, num_neg/(T - burn_in + 1))
		prob_0 = append(prob_0, num_0/(T - burn_in + 1))
	  
		if(num_0 > (T - burn_in + 1)/4){
			next;
		}
	    if(num_pos > num_neg){
	    	rho_mean[v] = 1
	    }
	    else if(num_pos < num_neg){
	    	rho_mean[v] = -1
	    }
	    else{
	    	next;
	    }
	}
	for(v in c(1:nrow(grids))){
		if(rho_true[v] > 0){
			if(rho_mean[v] > 0){
				sensi_pos_num = sensi_pos_num + 1
				speci_neg_num = speci_neg_num + 1
			}
			if(rho_mean[v] == 0){
				speci_neg_num = speci_neg_num + 1
			}
			# if(rho_mean[v] < 0){
			# 	FDR_neg_num = FDR_neg_num + 1
			# }
		}
		if(rho_true[v] == 0){
			if(rho_mean[v] == 0){
				speci_pos_num = speci_pos_num + 1
				speci_neg_num = speci_neg_num + 1
			}
			if(rho_mean[v] > 0){
				speci_neg_num = speci_neg_num + 1
				FDR_pos_num = FDR_pos_num + 1
			}
			if(rho_mean[v] < 0){
				speci_pos_num = speci_pos_num + 1
				FDR_neg_num = FDR_neg_num + 1
			}
			
		}
		if(rho_true[v] < 0){
			if(rho_mean[v] < 0){
				sensi_neg_num = sensi_neg_num + 1
				speci_pos_num = speci_pos_num + 1
			}
			if(rho_mean[v] == 0){
				speci_pos_num = speci_pos_num + 1
			}
			# if(rho_mean[v] > 0){
			# 	FDR_pos_num = FDR_pos_num + 1
			# }
		}
	}
	sensi_pos = sensi_pos_num/sum(rho_true == 1)
	sensi_neg = sensi_neg_num/sum(rho_true == -1)
	speci_pos = speci_pos_num/(nrow(grids) - sum(rho_true == 1))
	speci_neg = speci_neg_num/(nrow(grids) - sum(rho_true == -1))
	FDR_pos = FDR_pos_num /(sensi_pos_num + FDR_pos_num)
	FDR_neg = FDR_neg_num /(sensi_neg_num + FDR_neg_num)

	tau_1_sq_hat = rowSums(gibbs_tau_1_sq[,burn_in:T])/(T-burn_in +1)
	tau_2_sq_hat = rowSums(gibbs_tau_2_sq[,burn_in:T])/(T-burn_in +1)

	rho2 = rowSums(rho_hat)/(T-burn_in+1)
	xi_hat = rowSums(xi_hat_m[, c(burn_in:T)])/(T-burn_in +1)
	sigma_pos_sq_hat = rowSums(sigma_pos_sq_hat_m[, c(burn_in:T)])/(T-burn_in +1)
	sigma_neg_sq_hat = rowSums(sigma_neg_sq_hat_m[, c(burn_in:T)])/(T-burn_in +1)

	thres_xi = sigma_pos_sq_hat - sigma_neg_sq_hat
	abs_thres_xi = abs(thres_xi)
	var_Y_1_hat = tau_1_sq_hat + abs_thres_xi
	var_Y_2_hat = tau_2_sq_hat + abs_thres_xi
	cov_hat = sigma_pos_sq_hat - sigma_neg_sq_hat
	true_cov = dat$sigma_pos_sq - dat$sigma_neg_sq


	true_thres_xi = dat$sigma_pos_sq - dat$sigma_neg_sq
	true_abs_thres_xi = abs(true_thres_xi)
	true_var_Y_1 = dat$tau_1_sq + true_abs_thres_xi
	true_var_Y_2 = dat$tau_2_sq + true_abs_thres_xi
	fig_rho = twofigs.levelplot(
	 	rho_mean, dat$rho,
		grids[, 1],
	    grids[, 2],
	    titles = c("rho_hat", "rho")
	)
	fig_rho2 = twofigs.levelplot(
	 	rho2, dat$rho,
		grids[, 1],
	    grids[, 2],
	    titles = c("rho_hat", "rho")
	)
	fig_tau1 = twofigs.levelplot(
	 	tau_1_sq_hat, dat$tau_1_sq, 
		grids[, 1],
	    grids[, 2],
	    titles = c("tau_1_sq_hat", "tau_1_sq")
	)
	fig_tau2 = twofigs.levelplot(
	 	tau_2_sq_hat, dat$tau_2_sq, 
		grids[, 1],
	    grids[, 2],
	    titles = c("tau_2_sq_hat", "tau_2_sq")
	)
	fig_sigma_pos_sq = twofigs.levelplot(
     	sigma_pos_sq_hat, c(dat$sigma_pos_sq),
     	grids[, 1],
     	grids[, 2],
     	titles = c("sigma_pos_sq_hat", "sigma_pos_sq")
 	)

	fig_sigma_neg_sq = twofigs.levelplot(
     	sigma_neg_sq_hat, c(dat$sigma_neg_sq),
     	grids[, 1],
     	grids[, 2],
     	titles = c("sigma_neg_sq_hat", "sigma_neg_sq")
 	)
	
 	fig_cov = twofigs.levelplot(
     	cov_hat, true_cov,
     	grids[, 1],
     	grids[, 2],
     	titles = c("cov_hat", "cov")
 	)
 	fig_prob_map = threefigs.levelplot(
     	prob_pos, prob_neg, prob_0, 
     	grids[, 1],
     	grids[, 2],
     	titles = c("prob_pos", "prob_neg", "prob_0")
 	)
	return(list('xi_hat_m' = xi_hat_m, 'cor_type' = rho_mean, 'rho_hat' = rho_hat, 'var_Y_1_hat' = var_Y_1_hat, 'var_Y_2_hat' = var_Y_2_hat, 'cov_hat' = cov_hat, 
		'fig_rho' = fig_rho, 'fig_rho2' = fig_rho2,'fig_sigma_pos_sq' = fig_sigma_pos_sq, 'fig_sigma_neg_sq' = fig_sigma_neg_sq, 
		'fig_cov' = fig_cov, 'sensi_pos' = sensi_pos, 'sensi_neg' = sensi_neg, 'speci_pos' = speci_pos, 'speci_neg' = speci_neg, 'FDR_pos' = FDR_pos, 'FDR_neg' = FDR_neg, 'fig_tau1' = fig_tau1,
		'fig_tau2' = fig_tau2, 'fig_prob_map' = fig_prob_map, "prob_pos" = prob_pos, "prob_neg" = prob_neg, "prob_0" = prob_0))

}







get.chain.e = function(i, choice = 'pos'){
	if(choice == 'pos'){
		ch = matrix(nrow = L, ncol = T)
		for(t in c(1:T)){
			ch[,t] = gibbs_e_pos[[t]][,i]
		}
	}
	if(choice == 'neg'){
		ch = matrix(nrow = L, ncol = T)
		for(t in c(1:T)){
			ch[,t] = gibbs_e_neg[[t]][,i]
		}
	}
	return(t(ch))
}

