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
  real<lower = 0> alpha;
  real<lower = 0> r0;
  real ba1;
}

transformed parameters {
  real<lower=0, upper=1> prob[N];

  for (i in 1:N) {
    prob[i] = (1 - (1 - gamma_cdf(Ts1[i], alpha, r0)) ^ exp(ba1 * As[i]));
    prob[i] -= (1 - (1 - gamma_cdf(Ts0[i], alpha, r0)) ^ exp(ba1 * As[i]));
  }
}

model {
  r0 ~ exponential(1);
  ba1 ~ normal(0, 1);
  
  // target += sum(lp);
  target += log_sum_exp(prob);
  
  for (i in 1:N_Cen) {
    target += gamma_lccdf(Ts_Cen[i] | alpha, r0) * exp(ba1 * As_Cen[i]);
  }
}
