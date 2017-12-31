data {
  int<lower=1> N;                   // number of mutations

                                    // Signature-specific Data:
  int<lower=1> K;                   // number of reference signatures
  int<lower=1> R;                   // number of mutation types (vocabulary size, usually 96)
  int v[N];                         // mutation type dummy variable for each mutation
  matrix[K, R] ref_signatures;      // reference mutation signatures (S matrix in the model)

                                    // Clonality-specific Data:
  int<lower=1> L;                   // number of clones
  int x[N];                         // number of alt reads for each mutation
  int d[N];                         // number of total reads for each mutation (depth)
  vector[N] a;                      // correction factor (purity / (purity * Ct) + ((1 - purity) * Cn))
}

parameters {
  simplex[L * K] phi;         // mixing proportions - one coefficient per signature per clone

  ordered[L] mu_ordered; // subpopulation beta-binomial means
  real<lower=0> kappa_minus_two;    // beta-binomial concentration term
}

transformed parameters {
  vector<lower=0, upper=1>[L] mu = inv_logit(mu_ordered);
}

model {
  real alpha[N,L];       // beta-binomial: alpha = mu * kappa
  real beta[N,L];        // beta-binomial: beta = (1 - mu) * kappa
                            // Note that mu values are first adjusted by a values for each mutation

  int idx;                          // mapping l and k indices to corresponding coefficient of phi
  real sig_likelihood[N,K];         // Matrix of signature likelihoods for each mutation
  real pop_likelihood[N,L];         // Matrix of population likelihoods for each mutation
  real likelihood_sums[L, K, N];    // 

  kappa_minus_two ~ gamma(0.01, 0.01);
  mu_ordered ~ logistic(0, 1); // This is used instead of mu ~ beta(1, 1);

  alpha = to_array_2d((a * to_row_vector(mu)) * (kappa_minus_two + 2));
  beta = to_array_2d((1 - (a * to_row_vector(mu))) * (kappa_minus_two + 2));

  for (n in 1:N) {
    for (l in 1:L) {
      pop_likelihood[n,l] = beta_binomial_lpmf(x[n] | d[n], alpha[n,l], beta[n,l]);
    }
  }

  for (n in 1:N) {
    for (k in 1:K) {
      sig_likelihood[n,k] = categorical_lpmf(v[n] | to_vector(ref_signatures[k]));
    }
  }

  for (l in 1:L) {
    for (k in 1:K) {
      idx = ((l-1) * K) + k;
      likelihood_sums[l, k, :] = to_array_1d(log(phi[idx]) + to_vector(pop_likelihood[:,l]) + to_vector(sig_likelihood[:,k]));
    }
  }

  likelihood_sums = exp(likelihood_sums);

  for (n in 1:N) {
    target += log(sum(to_array_1d(likelihood_sums[:, :, n])));
  }
}

