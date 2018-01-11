data {
  int<lower=1> N; // number of mutations
  int<lower=1> L; // number of clones
  int x[N]; // number of alt reads for each mutation
  int d[N]; // number of total reads for each mutation
  vector<lower=1>[N] Ct; // tumour copy number at mutation site
  vector<lower=1>[N] Cn; // normal copy number at mutation site
  real purity; // tumour content as a fraction from 0 to 1
}
transformed data {
  vector[N] a = purity ./ ( (purity * Ct) + ((1 - purity) * Cn) );
}
parameters {
  simplex[L] clone_prop; // mixing proportions
  ordered[L] mu_ordered; // subpopulation beta-binomial means
  real<lower=0> kappa_minus_two;
}
transformed parameters {
  vector<lower=0, upper=1>[L] mu = inv_logit(mu_ordered);
}
model {
  real alpha[N,L];       // beta-binomial: alpha = mu * kappa
  real beta[N,L];        // beta-binomial: beta = (1 - mu) * kappa
                            // Note that mu values are first adjusted by a values for each mutation

  matrix[L,N] pop_likelihood;
  matrix[L,N] likelihood_sums;
  row_vector[N] likelihood_colsums;

  clone_prop ~ dirichlet(rep_vector(1, L));
  kappa_minus_two ~ gamma(0.01, 0.01);
  mu_ordered ~ logistic(0, 1); // This is used instead of mu ~ beta(1, 1);

  alpha = to_array_2d((a * to_row_vector(mu)) * (kappa_minus_two + 2));
  beta = to_array_2d((1 - (a * to_row_vector(mu))) * (kappa_minus_two + 2));

  for (n in 1:N) {
    for (l in 1:L) {
      pop_likelihood[l,n] = beta_binomial_lpmf(x[n] | d[n], alpha[n, l], beta[n, l]);
    }
  }

  for (l in 1:L) {
    likelihood_sums[l] = log(clone_prop[l]) + pop_likelihood[l];
  }

  likelihood_sums = exp(likelihood_sums);
  likelihood_colsums = to_row_vector(rep_vector(1, L)) * likelihood_sums;
  target += sum(log(likelihood_colsums));
}

