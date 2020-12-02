data {
  int<lower=1> G;                                                
  int<lower=1> S;                                                
  int<lower=0> expression[G, S];
  vector<lower=0, upper=1>[S] design;
  vector<lower=0>[G] log_l_g;
}

transformed data {
  vector<lower=-1, upper=1>[S] design2 = 2 * design - rep_vector(1, S);

}

parameters {
  real mu_theta;                       
  real<lower=0>  sigma_theta;          
  vector[G] theta;                     
  vector[G] beta;
  real<lower=0>  sigma_beta;
  vector[S] log_norm_factors;
  matrix[G, S] error;
  real<lower=0> sigma_error;
}

model {
  mu_theta ~ normal(5, 1.5) ;
  sigma_theta ~ std_normal();
  theta ~ normal(mu_theta, sigma_theta);
  sigma_beta ~ std_normal();
  beta ~ normal(0, sigma_beta);
  to_row_vector(error) ~ normal(0, sigma_error);
  sigma_error ~ normal(0, 1);
  log_norm_factors ~ normal(0, 0.05);
  
  for (i in 1:G) {
    for (j in 1:S) {
      expression[i, j] ~ poisson_log(theta[i] + beta[i] * design2[j] + error[i, j] + log_l_g[i] + log_norm_factors[j]);
    }
  }
}
