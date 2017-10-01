#' t-SNE plot of MCMC density clusters
#'
#' t-SNE is a useful tool for visualizing data clusters in high dimensional space.
#' This method produces a t-SNE plot showing the relationship of points sampled
#' by MCMC and coloured by density clusters.
#'
#' @param exposures_mcmc_output         Output from get_exposures
#' @param patient_density_clustering    Output from density_clustering
#'
#' @return A ggplot object containing the t-SNE plot in 2 dimensions
#'
#' @import ggplot
#' @import dplyr
#' @import tibble
#' @import Rtsne
#' 
#' @export

plot_density_tsne <- function(exposures_mcmc_output, patient_density_clustered, perplexity = 30) {
  mcmc_chain <- exposures_mcmc_output %>% 
    merge_exposures_with_clustering(patient_density_clustered) %>%
    group_by_at(exposures_mcmc_output$signature_names) %>%
    filter(row_number() == 1) %>%
    ungroup()
  
  tsne_data <- mcmc_chain %>%
    select(-iteration, -chain, -cluster) %>%
    as.matrix() %>%
    Rtsne(
        perplexity = perplexity,
        max_iter = 5000
    )
  
  mcmc_chain %>%
    mutate(
      tsne_1 = tsne_data$Y[, 1],
      tsne_2 = tsne_data$Y[, 2],
      cluster = cluster %>% as.factor %>% 
        recode(
          `0` = 'Outlier'
        )
    ) %>% 
    select(chain, iteration, cluster, tsne_1, tsne_2) %>%
    ggplot(aes(
      x = tsne_1,
      y = tsne_2,
      colour = cluster,
      group = chain
    )) +
    geom_path(colour = 'grey80') +
    geom_point()
}
