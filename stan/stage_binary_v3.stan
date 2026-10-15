data {
  int<lower=1> N;
  int<lower=0> K;
  matrix[N, K] X;
  int<lower=1> G;
  array[G] int<lower=1> n_level;
  array[N, G] int<lower=1> idx;
  array[N] int<lower=1> n;
  array[N] int<lower=0> y;
  int<lower=0, upper=2> use_t;   // 0: 正規, 1: t（自由度を推定）, 2: t（自由度を固定）
  real<lower=1> nu_fixed;
}
transformed data {
  int n_total = sum(n_level);
  array[G, N] int pos;
  {
    int s = 0;
    for (g in 1:G) {
      for (i in 1:N) pos[g, i] = s + idx[i, g];
      s += n_level[g];
    }
  }
}
parameters {
  real alpha;
  vector[K] beta;
  vector<lower=0>[G] sigma;
  vector[n_total] z;
  array[use_t == 1 ? 1 : 0] real<lower=1> nu;
}
model {
  vector[N] eta = rep_vector(alpha, N);
  if (K > 0) eta += X * beta;
  for (g in 1:G) eta += sigma[g] * z[pos[g]];

  alpha ~ normal(0, 1.5);
  beta ~ normal(0, 1);
  sigma ~ normal(0, 1);
  if (use_t == 0) {
    z ~ std_normal();
  } else {
    if (use_t == 1) {
      nu[1] ~ gamma(2, 0.1);
      z[1:n_level[1]] ~ student_t(nu[1], 0, 1);
    } else {
      z[1:n_level[1]] ~ student_t(nu_fixed, 0, 1);
    }
    z[(n_level[1] + 1):n_total] ~ std_normal();
  }
  y ~ binomial_logit(n, eta);
}
