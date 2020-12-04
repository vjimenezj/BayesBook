data {
  int<lower=1> G;                     
  int<lower=1> S;                     
  vector<lower=0, upper=1>[S] design;
}

transformed data {
  vector<lower=-1, upper=1>[S] design2 = 2 * design - rep_vector(1, S);
}

generated quantities {
  real  mu_alpha = 5;                 
  real<lower=0>  sigma_alpha = 1.5;     
  real<lower=0>  sigma_beta = 1;     
  vector[G] alpha;                    
  vector[G] beta;                    
  int<lower=0> expression[G, S]; 
  matrix[G, S] error;
  for (i in 1:G) {
    alpha[i] = normal_rng(mu_alpha, sigma_alpha);
  }
  for (i in 1:G) {
    beta[i] = normal_rng(0, sigma_beta);
  }
  for (i in 1:G) {
    for (j in 1:S) {
    error[i, j] = normal_rng(0, 0.3);
    expression[i, j] = poisson_log_rng(alpha[i] + beta[i] * design2[j] + error[i, j]);
    }
  }
}
