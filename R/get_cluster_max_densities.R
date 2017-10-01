#' Reports the maximum density point in each cluster
#'
#' As MCMC chains converge, sampling density becomes proportional to the posterior distribution.
#' Thus, areas of greater density in the MCMC solution space may represent more promising solutions.
#' As a simple heuristic for density, this function returns the point with the least mean distance
#' to its nearest k neighbours. By default \eqn{k=400}, but the user may wish to adjust this value
#' to be lower if the number of iterations in decreased, or higher if the number of iterations is
#' increased.
#'
#' @param exposures_mcmc_output     Output from plot_exposures()
#' @param exposure_mcmc_clusters    Output from density_clustering()
#' @param k_value                   The value of k (default 400)
#'
#' @return A dataframe of solutions, one per cluster, in long format (one row per signature per cluster)
#'
#' @import dplyr
#' @import tibble
#' @import tidyr

get_cluster_max_densities <- function(exposures_mcmc_output, exposure_mcmc_clusters, k_value = 400) {
  exposures_matrix <- mcmc_exposures_as_matrix(exposures_mcmc_output)
  
  max_density_per_cluster <- merge_exposures_with_clustering(exposures_mcmc_output, exposure_mcmc_clusters) %>%
    mutate(mean_knn_dist = kNNdist(exposures_matrix, k = k_value) %>% apply(1, mean)) %>%
    group_by(cluster) %>%
    filter(mean_knn_dist == max(mean_knn_dist)) %>%
    ungroup() %>%
    left_join(
      tibble(
        cluster_score = exposure_mcmc_clusters$cluster_scores
      ) %>%
        mutate(cluster = row_number()),
      by = c('cluster')
    ) %>%
    gather(signature, exposure, -chain, -iteration, -cluster, -mean_knn_dist, -cluster_score)
  
  return(max_density_per_cluster)
}
