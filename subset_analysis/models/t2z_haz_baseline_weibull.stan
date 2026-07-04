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
  real<lower=0> alpha;
  real<lower=0> sigma;
  real ba0;
  real ba1;
}

transformed parameters {
  real<lower=0, upper=1> prob[N];

  for (i in 1:N) {
    prob[i] = (1 - (1 - weibull_cdf(Ts1[i], alpha, sigma)) ^ exp(ba0 + ba1 * As[i]));
    prob[i] -= (1 - (1 - weibull_cdf(Ts0[i], alpha, sigma)) ^ exp(ba0 + ba1 * As[i]));
  }
}

model {
  alpha ~ exponential(1);
  sigma ~ exponential(1);
  ba0 ~ normal(0, 1);
  ba1 ~ normal(0, 1);
  
  // target += sum(lp);
  target += log_sum_exp(prob);
  
  for (i in 1:N_Cen) {
    target += weibull_lccdf(Ts_Cen[i] | alpha, sigma) * exp(ba0 + ba1 * As_Cen[i]);
  }
}

generated quantities {
  real<lower = 0> tte_sim;
  
  tte_sim = weibull_rng(alpha, sigma);
}
