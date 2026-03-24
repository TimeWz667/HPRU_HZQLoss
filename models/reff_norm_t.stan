
data {
    int N;
    int N_gp; // number of groups

    real Ys[N]; // outcome
    real Ts[N]; // outcome
    int Gs[N]; // group 
}
parameters {
    real b0;
    real bt;
    real bg[N_gp]; 

    real<lower=0,upper=10> sigma;  
}
model {
  bg ~ normal(0, 0.1);
  b0 ~ normal(0, 1);
  bt ~ normal(0, 1);

  for(i in 1:N) {
    target += normal_lpdf(Ys[i] | b0 + bt * Ts[i] + bg[Gs[i]], sigma);
  }
}
