data {
  int<lower=1> G;                                                
  int<lower=1> S;                                                
  int<lower=0> expression[G, S];
  vector<lower=0, upper=1>[S] design;
  vector<lower=0>[G] eff_length;
}

transformed data {
  vector<lower=-1, upper=1>[S] design2 = 2 * design - rep_vector(1, S);       
  real <lower=0> expected_gene_mean = sum(expression[, 1])/sum(eff_length);
  vector<lower=0>[G] log_eff_length = log(eff_length);
}

parameters {
  real mu_alpha;                       
  real<lower=0>  sigma_alpha;          
  vector[G] alpha;                     
  vector[G] beta;
  real<lower=0>  sigma_beta;
  vector[S-1] log_norm_factors_ref;
  real<lower=0, upper=2000000000> phi;
}

transformed parameters {
  vector[S] log_norm_factors = append_row(0, log_norm_factors_ref);
}

model {
  mu_alpha ~ normal(log(expected_gene_mean), 1.17) ;
  sigma_alpha ~ normal(0, 2);
  alpha ~ normal(mu_alpha, sigma_alpha);
  sigma_beta ~ std_normal();
  beta ~ normal(0, sigma_beta);
  log_norm_factors_ref ~ normal(0, 0.05);
  phi ~ uniform(0, 2000000000);
  for (i in 1:G) {
    for (j in 1:S) {
      expression[i, j] ~ neg_binomial_2_log(alpha[i] + beta[i] * design2[j] + log_eff_length[i] + log_norm_factors[j], phi);
    }
  }
}

generated quantities {
  matrix[G, S] log_lik;                                                           
  int expression_rep[G, S]; 
  for (i in 1:G) {
    for (j in 1:S) {
      log_lik[i, j] = neg_binomial_2_log_lpmf(expression[i, j] | alpha[i] + beta[i] * design[j] + log_norm_factors[j] + log_eff_length[i], phi);
      expression_rep[i, j] = neg_binomial_2_log_rng(alpha[i] + beta[i] * design[j] + log_norm_factors[j] + log_eff_length[i], phi);
    }
  }
}
