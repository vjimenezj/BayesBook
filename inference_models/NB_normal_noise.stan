data {
  int N;
  vector[N] alpha;
  real<lower=0> phi;
  real<lower=0> sigma_error;
}

parameters {
  vector<lower=0>[N] gamma;
  vector<lower=0>[N] error;
  real<lower=0, upper=2000000000> mixed_phi;
}

transformed parameters{
  vector<lower=0>[N] mixed;
  mixed = gamma + error;
}  

model{
  for (i in 1:N) {
    gamma[i] ~ gamma(phi, phi/alpha[i]) ;
  }
  for (i in 1:N) {
    error[i] ~ lognormal(0, sigma_error);
  }
  mixed_phi ~ uniform(0, 2000000000);
  for (i in 1:N) {
    mixed[i] ~ gamma(mixed_phi, mixed_phi/alpha[i]);
  }
}

