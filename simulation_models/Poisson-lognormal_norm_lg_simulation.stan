data {
  int<lower=1> G;                     
  int<lower=1> S;                     
  vector<lower=0, upper=1>[S] design;
  vector<lower=0>[G] log_l_g;
}

transformed data {
  vector<lower=-1, upper=1>[S] design2 = 2 * design - rep_vector(1, S);
}

generated quantities {
  real  mu_alpha = 6;                 
  real<lower=0>  sigma_alpha = 1.5;     
  real<lower=0>  sigma_beta = 1;     
  real<lower=0>  sigma_error = 0.5;     
  vector[G] alpha;                    
  vector[G] beta;                    
  real expression[G, S]; 
  vector[S] log_norm_factors;
  matrix[G, S] error;
  for (i in 1:G) {
    alpha[i] = normal_rng(mu_alpha, sigma_alpha);
  }
  log_norm_factors[1] = 0;
  for (i in 2:S) {
    log_norm_factors[i] = normal_rng(0, 0.05);
  }
  for (i in 1:G) {
    beta[i] = normal_rng(0, sigma_beta);
  }
  for (i in 1:G) {
    for (j in 1:S) {
    error[i, j] = normal_rng(0, sigma_error);
    expression[i, j] = poisson_log_rng(alpha[i] + beta[i] * design2[j] + error[i, j] + log_norm_factors[j] + log_l_g[i]);
    }
  }
}
