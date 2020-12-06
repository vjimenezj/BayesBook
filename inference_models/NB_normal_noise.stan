data {
  real<lower=0> phi;
  real<lower=0> sigma_error;
}

transformed data {
  int N=1000;
  int M=10000;
}

parameters {
  vector[N] alpha;
  matrix[N, M] gamma;
  matrix[N, M] error;
  real<lower=0, upper=2000000000> mixed_phi;
}

transformed parameters{
  matrix[N, M] mixed;
  mixed = gamma + error;
}  

model{
  alpha ~ normal(5, 1.5);
  for (i in 1:N) {
    gamma[i, ] ~ gamma(phi, phi/alpha[i]) ;
  }
  for (i in 1:N) {
    error[i, ] ~ normal(0, sigma_error);
  }
  mixed_phi ~ uniform(0, 2000000000);
  for (i in 1:N) {
    mixed[i, ] ~ gamma(mixed_phi, mixed_phi/alpha[i]);
  }
}

