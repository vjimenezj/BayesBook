data {
  int N;
  vector[N] alpha;
  real<lower=0> phi;
  real<lower=0> sigma_error;
}

transformed data {
  int M=1000;
}

parameters {
  matrix[N, M] gamma;
  matrix[N, M] error;
  real<lower=0, upper=2000000000> mixed_phi;
}

transformed parameters{
  matrix[N, M] mixed;
  mixed = gamma + error;
}  

model{
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

