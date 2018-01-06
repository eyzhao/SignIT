#' Plots a simple NNLS solution
#'
#' Obtains the naive NNLS signature solution and plots the exposures.
#'
#' Non-negative least squares computes the value of x which minimizes the residual sum
#' of squares \eqn{\lVert Ax = b \rVert^2} given a known matrix A and vector b.
#' This provides a simple solution to the mutation signature problem. However,
#' merely minimizing the residual sum of squares may result in overfitting the solution.
#' This plotting function computes and plots a simple NNLS-based solution to the
#' exposure vector and returns a plot. This can serve as a baseline comparison or sanity check
#' against the MCMC methods featured in SignIT.
#'
#' @param exposures_mcmc_output
#'
#' @return A ggplot object showing the NNLS solution for each exposure value.
#'
#' @import nnls
#' @import tibble
#' @import dplyr
#' 
#' @export

plot_nnls_solution <- function(exposures_mcmc_output, signature_trim='Signature') {
  nnls_solution <- nnls(
    exposures_mcmc_output$reference_signatures %>% select(-mutation_type) %>% as.matrix,
    exposures_mcmc_output$mutation_catalog$count
  )$x
  
  nnls_plot <- tibble(
      signature = exposures_mcmc_output$signature_names %>% factor(levels = exposures_mcmc_output$signature_names),
      nnls = nnls_solution
    ) %>%
    mutate(signature = trim_signature_names(signature, signature_trim)) %>%
    ggplot(aes(x = signature, y = nnls)) +
    geom_point() +
    rotate_x_axis_labels()
}


