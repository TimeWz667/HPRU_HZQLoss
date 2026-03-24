data {
  int<lower=0> N;
  real Ts[N];
  int As[N];
}

parameters {
  real<lower = 0> r0;
  real ba1;
}

model {
  r0 ~ exponential(1);
  ba1 ~ normal(0, 1);

  for (i in 1:N) {
    target += exponential_lpdf(Ts[i] | r0 * exp(ba1 * As[i]));
  }
}
