#' Get the STAN Model for Mutation Signature Analysis
#'
#' Creates the STAN model for MCMC analysis of mutation signatures.
#'
#' @return STAN model object
#'
#' @import rstan
#' @import nnls
#' @import dplyr

get_signature_model <- function() {

    stan_model_string_vectorized <- '
        data {
          int<lower=1> S; // number of signatures
          int<lower=1> N; // number of mutations
          int<lower=1> R; // number of mutation types (vocabulary size)
          int<lower=0> v[R]; // mutation catalog vector
          matrix[R, S] ref_signatures; // reference mutation signatures (theta)
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
          v ~ multinomial(sim_catalog_prob);
          exposures ~ dirichlet(rep_vector(1, S));
        }

    '

    message('Establishing Stan model')

    stan_dso = stan_model(model_code = stan_model_string_vectorized)

}
