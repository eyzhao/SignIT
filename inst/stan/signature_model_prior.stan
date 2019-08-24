data {
  int<lower=1> S; // number of signatures
  int<lower=1> N; // number of mutations
  int<lower=1> R; // number of mutation types (vocabulary size)
  int<lower=0> v[R]; // mutation catalog vector
  matrix[R, S] ref_signatures; // reference mutation signatures (theta)
  real<lower=0,upper=1> prior[S];
}

parameters {
  simplex[S] exposures; // mixing proportions
}

transformed parameters {
  vector[R] sim_catalog;
  simplex[R] sim_catalog_prob;
  sim_catalog = ref_signatures * exposures * N;
  sim_catalog_prob = sim_catalog / sum(sim_catalog);
}

model {
  for (s in 1:S) {
    if (exposures[s] < 0.001) {
      target += log((1 - prior[s])/0.001);
    } else {
      target += log(prior[s] / 0.999);
    }
  }
  v ~ multinomial(sim_catalog_prob);
}

generated quantities {
  real log_lik;

  log_lik = multinomial_lpmf(v | sim_catalog_prob);
}

