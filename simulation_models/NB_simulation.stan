data {
  int<lower=1> G;                                                
  int<lower=1> S;                                                
  vector<lower=0, upper=1>[S] design;
  vector<lower=0>[G] eff_length;
  real<lower=0, upper=2000000000> phi;
}

transformed data {
  vector<lower=-1, upper=1>[S] design2 = 2 * design - rep_vector(1, S);       
  vector<lower=0>[G] log_eff_length = log(eff_length);
}

generated quantities {
  real  mu_alpha = 2;                 
  real<lower=0>  sigma_alpha = 2;     
  real<lower=0>  sigma_beta = 1;     
  vector[G] alpha;                    
  vector[G] beta;                    
  int<lower=0> expression[G, S]; 
  vector[S] log_norm_factors;
  for (i in 1:G) {
    alpha[i] = normal_rng(mu_alpha, sigma_alpha);
  }
  for (i in 1:G) {
    beta[i] = normal_rng(0, sigma_beta);
  }
  log_norm_factors[1] = 0;
  for (i in 2:S) {
    log_norm_factors[i] = normal_rng(0, 0.05);
  }
  for (i in 1:G) {
    for (j in 1:S) {
      expression[i, j] = neg_binomial_2_log_rng(alpha[i] + beta[i] * design2[j] + log_eff_length[i] + log_norm_factors[j], phi);
    }
  }
}
