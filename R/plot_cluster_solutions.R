#' Plots solutions derived from density clustering
#'
#' This plots the output of get_cluster_summary_points. See get_cluster_summary_points.R for
#' more details on the rationale behind this plot.
#'
#' @param exposures_mcmc_output     Output from get_exposures
#' @param exposure_mcmc_clusters    Output from density_clustering
#' @param method                    Either 'bestfit' or 'mostdense'

plot_cluster_solutions <- function(exposures_mcmc_output, exposure_mcmc_clusters, method = 'bestfit', k_value = 400) {
  cluster_summary_points <- get_cluster_summary_points(
    exposures_mcmc_output, 
    exposure_mcmc_clusters, 
    method, 
    k_value
  ) %>%
    mutate(signature = signature %>% factor(levels = exposures_mcmc_output$signature_names))
  
  cluster_summary_points %>%
    mutate(signature = trim_signature_names(signature)) %>%
    ggplot(aes(
      x = signature,
      y = exposure
    )) +
    geom_line(aes(
      group = cluster,
      alpha = cluster_score
    )) +
    rotate_x_axis_labels()
    
}


