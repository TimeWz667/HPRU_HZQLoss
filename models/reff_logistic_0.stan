
data {
    int N;
    int N_gp; // number of groups

    int Ys[N]; // outcome
    int Gs[N]; // group 
}
parameters {
    real b0;
    real bg[N_gp]; 
}
model {
  bg ~ normal(0, 0.1);
  b0 ~ normal(0, 5);

  for(i in 1:N) {
    target += bernoulli_lpmf(Ys[i] | inv_logit(b0 + bg[Gs[i]]));
  }
}
