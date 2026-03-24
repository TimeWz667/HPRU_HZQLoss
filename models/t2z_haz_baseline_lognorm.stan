data {
  int<lower=0> N;
  real Ts0[N];
  real Ts1[N];
  int As[N];
  
  int<lower=0> N_Cen;
  real Ts_Cen[N_Cen];
  int As_Cen[N_Cen];
}

parameters {
  // hyper parameters
  real mu;
  real<lower=0> sigma;
  real ba0;
  real ba1;
}

transformed parameters {
  real<lower=0, upper=1> prob[N];
  real<lower=0> mu_e;
  mu_e = exp(mu + sigma^2/2);

  for (i in 1:N) {
    prob[i] = (1 - (1 - lognormal_cdf(Ts1[i], mu, sigma)) ^ exp(ba0 + ba1 * As[i]));
    prob[i] -= (1 - (1 - lognormal_cdf(Ts0[i], mu, sigma)) ^ exp(ba0 + ba1 * As[i]));
  }
}

model {
  ba1 ~ normal(0, 1);
  
  target += normal_lpdf(mu_e | 10, 2);
  
  // target += sum(lp);
  target += log_sum_exp(prob);
  
  for (i in 1:N_Cen) {
    target += lognormal_lccdf(Ts_Cen[i] | mu, sigma) * exp(ba0 + ba1 * As_Cen[i]);
  }
}

generated quantities {
  real<lower = 0> tte_sim;
  
  tte_sim = lognormal_rng(mu, sigma);
}
