#' Get the STAN Model for Mutation Signature Analysis
#'
#' Creates the STAN model for MCMC analysis of mutation signatures.
#'
#' @return STAN model object
#'
#' @import rstan
#' @import nnls
#' @import dplyr

get_stan_model <- function(model_type) {
    message('Establishing Stan model')

    if (model_type == 'signature') {
        model_path = paste0(path.package('signit'), '/inst/signature_model.stan')
    } else if (model_type == 'population') {
        model_path = paste0(path.package('signit'), '/inst/population_model.stan')
    } else if (model_type == 'joint') {
        model_path = paste0(path.package('signit'), '/inst/signit_model_infer_populations.stan')
    }

    stan_dso = stan_model(model_path)
}
