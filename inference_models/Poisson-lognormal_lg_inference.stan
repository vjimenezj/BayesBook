data {
  int<lower=1> G;                                                
  int<lower=1> S;                                                
  int<lower=0> expression[G, S];
  vector<lower=0, upper=1>[S] design;
  vector<lower=0>[G] log_l_g;
}

transformed data {
  vector<lower=-1, upper=1>[S] design2 = 2 * design - rep_vector(1, S);
  real<lower=0> mean_eff_length = mean(log_l_g);
  vector[G] eff_log_l_g = log_l_g - rep_vector(mean_eff_length, G);
  real<lower=0> expected_mean;
  for (i in 1:S) {
    expected_mean = sum(expression[, i]) * 1.0 / (G * S);
  }
}

parameters {
  real mu_alpha;                       
  real<lower=0>  sigma_alpha;          
  vector[G] alpha;                     
  vector[G] beta;
  real<lower=0>  sigma_beta;
  matrix[G, S] error;
  real<lower=0> sigma_error;
}

model {
  mu_alpha ~ normal(expected_mean, 1.5) ;
  sigma_alpha ~ std_normal();
  alpha ~ normal(mu_alpha, sigma_alpha);
  sigma_beta ~ std_normal();
  beta ~ normal(0, sigma_beta);
  to_row_vector(error) ~ normal(0, sigma_error);
  sigma_error ~ std_normal();

  for (i in 1:G) {
    for (j in 1:S) {
      expression[i, j] ~ poisson_log(alpha[i] + beta[i] * design2[j] + error[i, j] + eff_log_l_g[i]);
    }
  }
}
