
data {
    int N;
    int N_gp; // number of groups

    int Ys[N]; // outcome
    real Ts[N]; // outcome
    real As[N];
    int Gs[N]; // group 
}
parameters {
    real b0;
    real bd15;
    real bd30;
    real ba1;
    real ba2;
    real bg[N_gp];
}

transformed parameters {
  real d15[N];
  real d30[N];
  real as[N];
  
  for (i in 1:N) {
    d15[i] = (Ts[i] < 15 / 365.25)?1:0;
    d30[i] = (Ts[i] < 30 / 365.25)?1:0;
    as[i] = (As[i] - 75) / 50.0;
  }
}

model {
  bg ~ normal(0, 0.1);
  b0 ~ normal(0, 5);
  bd15 ~ normal(0, 1);
  bd30 ~ normal(0, 1);
  ba1 ~ normal(0, 1);
  ba2 ~ normal(0, 1);

  for(i in 1:N) {
    target += bernoulli_lpmf(Ys[i] | inv_logit(b0 + d15[i] * bd15 + d30[i] * bd30 + as[i] * ba1 + (as[i] ^ 2) * ba2 + bg[Gs[i]]));
  }
}
