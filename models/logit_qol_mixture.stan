data {
  int<lower=0> N;
  real Ys[N];
}

transformed data {
  real logit_Ys[N];
  logit_Ys = logit(Ys);
}

parameters {
  real<lower=0, upper=1> p0;
  real mu[2];
  real<lower=0> sigma[2];
}

model {
  p0 ~ beta(1, 1);
  sigma ~ gamma(1, 1);
  
  target += log_sum_exp(log(p0) + normal_lpdf(Ys | mu[1], sigma[1]), log1m(p0) + normal_lpdf(Ys | mu[2], sigma[2]));
}
