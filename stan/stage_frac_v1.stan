data {
  int<lower=1> N;
  int<lower=0> K;
  matrix[N, K] X;
  int<lower=1> G;
  array[G] int<lower=1> n_level;
  array[N, G] int<lower=1> idx;
  vector<lower=0>[N] n;          // 打球数
  vector<lower=0>[N] ystar;      // 中立化した確率の合計
  real<lower=0> phi;             // 情報量の比（尤度の重み）
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
}
model {
  vector[N] eta = rep_vector(alpha, N);
  if (K > 0) eta += X * beta;
  for (g in 1:G) eta += sigma[g] * z[pos[g]];

  alpha ~ normal(0, 1.5);
  beta ~ normal(0, 1);
  sigma ~ normal(0, 1);
  z ~ std_normal();
  target += phi * (dot_product(ystar, log_inv_logit(eta)) +
                   dot_product(n - ystar, log1m_inv_logit(eta)));
}
