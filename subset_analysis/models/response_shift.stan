data {
  int N;
  int Ts[N];
  real Qs[N];
}
parameters {
  real qol0;
  real qol_diff;
  real<lower = 0> qol_drop;
  real<lower = 0> rec;
  real<lower = 0> err;
  
}
transformed parameters {
  real qol_final = qol0 - qol_diff;
  real qol_dz = qol0 - qol_drop;
  real red = qol_final - qol_dz;
}
model {
  for (i in 1:N) {
    
    if (Ts[i] < 0) {
      target += normal_lpdf(Qs[i]| qol0, err);
    } else {
      target += normal_lpdf(Qs[i]| qol_final - red * exp(- rec * (Ts[i] - 0)), err);
    }
   
  }
}
