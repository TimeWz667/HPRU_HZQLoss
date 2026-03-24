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
  real<lower = 0> r0;
  real ba1;
}

transformed parameters {
  // real<lower=0, upper=1> prob[N];
  // 
  // for (i in 1:N) {
  //   prob[i] = exponential_cdf(Ts1[i], r0 * exp(ba1 * As[i])) - exponential_cdf(Ts0[i], r0 * exp(ba1 * As[i]));
  // }
  real lp[N];
  real rate;
  
  for (i in 1:N) {
    rate = r0 * exp(ba1 * As[i]);
    lp[i] = - rate * Ts0[i] + exponential_lcdf(Ts1[i] - Ts0[i]| rate);
  }
  
  
}

model {
  r0 ~ exponential(1);
  ba1 ~ normal(0, 1);
  
  target += sum(lp);
  
  for (i in 1:N_Cen) {
    target += exponential_lccdf(Ts_Cen[i] | r0 * exp(ba1 * As_Cen[i]));
  }
}
