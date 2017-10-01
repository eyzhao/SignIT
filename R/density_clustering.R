#' Cluster MCMC solutions by density
#'
#' Uses the HDBSCAN method to perform density clustering. This analysis can be
#' used on the MCMC solution space to identify density clusters which may
#' represent alternative solutions of interest.
#'
#' @param exposures_mcmc_output     Output from get_exposures()
#'
#' @param minPts                    Minimum number of points in a neighbourhood
#'                                  for inclusion in a density cluster. The default
#'                                  is 40. Increasing this value will result in
#'                                  more clusters and more outliers.
#'
#' @param feature_selection_threshold   For each signature, the percentage of total mutation
#'                                      burden is computed across all MCMC iterations. If the
#'                                      maximum exposure is less than the feature selection
#'                                      threshold, then the signature is removed from density
#'                                      clustering.
#'
#' @return The HDBSCAN output, which contains data about the clusters corresponding to MCMC exposures
#'
#' @importFrom dbscan hdbscan
#'
#' @export

density_clustering <- function(exposures_mcmc_output, minPts = 50, feature_selection_threshold = 0.05) {
  suprathreshold_signatures <- exposures_mcmc_output$exposure_chain %>%
    mutate(exposure = exposure / exposures_mcmc_output$n_mutations) %>%
    group_by(signature) %>%
    summarise(
      max_exposure = max(exposure)
    ) %>%
    filter(max_exposure > feature_selection_threshold) %>%
    .$signature %>%
    as.character
  
  exposures_matrix <- mcmc_exposures_as_matrix(exposures_mcmc_output) 
  exposures_matrix <- exposures_matrix %>%
    as_tibble() %>%
    .[, suprathreshold_signatures] %>%
    as.matrix()
  
  # This will take a few minutes
  patient_density_clustered <- hdbscan(exposures_matrix, minPts=minPts)
  
  return(patient_density_clustered)
}
