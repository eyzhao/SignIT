#' Report the point in each cluster with the lowest residual sum of squares
#'
#' Given the known mutation catalog (m), a reference signature matrix (S) and
#' a potential exposure vector (e), the residual sum of squares (RSS) is defined as
#'
#' \deqn{\lVert \mathbf{m} - \mathbf{Se} \rVert^2}
#'
#' such that minimizing RSS produces the "best fit" in terms of Euclidean distance
#' between between the actual mutation catalog and the catalog simulated using the
#' computed exposure matrix.
#'
#' @param exposures_mcmc_output     Output from get_exposures() function.
#' @param exposure_mcmc_clusters    Output from density_clustering() function.
#'
#' @return A dataframe of solutions, one per cluster, in long format (one row per signature per cluster)
#'
#' @import dplyr
#' @import tibble
#' @import tidyr

get_cluster_best_fits <- function(exposures_mcmc_output, exposure_mcmc_clusters) {
  best_fit_per_cluster <- merge_exposures_with_clustering(exposures_mcmc_output, exposure_mcmc_clusters) %>%
    mutate(rss = get_rss(exposures_mcmc_output)) %>%
    group_by(cluster) %>%
    filter(rss == min(rss)) %>%
    ungroup() %>%
    select(-iteration, -chain) %>%
    distinct() %>%
    left_join(
      tibble(
        cluster_score = exposure_mcmc_clusters$cluster_scores
      ) %>%
        mutate(cluster = row_number()),
      by = c('cluster')
    ) %>%
    gather(signature, exposure, -cluster, -rss, -cluster_score)
  
  return(best_fit_per_cluster)
}


