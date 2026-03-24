data {
  int<lower=0> N;
  int Ys[N];
  real Ts[N];
}

parameters {
  real<lower=0> rate;
}

transformed parameters {
  real<lower=0, upper=1> prob[N];
  
  for (i in 1:N) {
    prob[i] = exponential_cdf(Ts[i], rate);
  }
  
}

model {
  rate ~ exponential(1);
  
  for (i in 1:N) {
    target += binomial_lpmf(Ys[i] | 1, prob[i]);
  }
}
