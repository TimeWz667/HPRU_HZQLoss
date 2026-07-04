// source: https://mc-stan.org/docs/2_18/stan-users-guide/soft-k-means.html
data {
  int<lower=0> N;
  int<lower=1> K;
  real Ys[N];
}
transformed data {
  real<upper=0> neg_log_K;
  neg_log_K = -log(K);
}
parameters {
  real mu[K]; // cluster means
}
transformed parameters {
  real<upper=0> soft_z[N, K]; // log unnormalized clusters
  for (n in 1:N) {
    for (k in 1:K) {
      soft_z[n, k] = neg_log_K - 0.5 * (mu[k] - Ys[n]) ^ 2;
    }
  }
}
model {
  // prior
  for (k in 1:K)
    mu[k] ~ std_normal();

  // likelihood
  for (n in 1:N)
    target += log_sum_exp(soft_z[n]);
}
