#' Summary plot of clusters
#'
#' @param exposures_mcmc_output     Output from plot_exposures
#' @param exposure_mcmc_clusters    Output from density_clustering
#' @param method                    Either 'bestfit' or 'mostdense'
#' 
#' @return ggplot object. Scatterplot of clusters showing sum of squared errors vs. cluster score.
#'
#' @import dplyr
#' @import tibble
#' @import ggplot
#'
#' @export

plot_cluster_solution_errors <- function(exposures_mcmc_output, exposure_mcmc_clusters) {
  cluster_score_table <- tibble(
    cluster = names(exposure_mcmc_clusters$cluster_scores) %>% as.numeric,
    cluster_score = exposure_mcmc_clusters$cluster_scores
  )
  
  cluster_summary_points <- merge_exposures_with_clustering(
    exposures_mcmc_output, 
    exposure_mcmc_clusters
  ) %>%
    mutate(rss = get_nnls_rss(exposures_mcmc_output)) %>%
    inner_join(cluster_score_table, by = 'cluster') %>%
    group_by(cluster) %>%
    summarise(
      cluster_score = mean(cluster_score),
      rss_mean = mean(rss)
    ) %>%
    ggplot(aes(
      x = cluster_score,
      y = rss_mean,
      ymin = rss_lCI,
      ymax = rss_uCI
    )) +
    geom_point()
}
