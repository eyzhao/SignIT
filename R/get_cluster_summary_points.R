#' Summarize each cluster using representative points
#'
#' The MCMC chains explore potential solutions in high dimensional space, where each signature
#' exposure value can vary. For the 30-signature reference, this amounts to a 30-dimensional space.
#' This is difficult to visualize. SignIT employs density clustering to explore this high dimensional space. 
#' However, it is similarly difficult to visualize the space of density clusters. 
#'
#' One approach is to select one representative point per density cluster which provides the best
#' solution according to some predetermined metric. This function reduces the exposure MCMC solution space
#' to one such solution per cluster, which can be plotted to explore a small subset of good alternative solutions.
#'
#' This function provides an interface to the various "summarization" methods. At present, two methods are available.
#' The first (default) method is "bestfit", which identifies the one point in each density cluster which minimizes
#' the residual sum of squares, and thus reduces the Euclidean distance between the solution and the input
#' mutation catalog. The second method is "mostdense", which estimates the highest density point within each cluster.
#'
#' This method is provided for a tabular representation of the summary points. For a graphical view,
#' we also provide the method plot_cluster_solutions().
#'
#' @param exposures_mcmc_output   Output from get_exposures()
#' 
#' @param exposure_mcmc_clusters  Output from density_clustering()
#'
#' @param method                  Either "bestfit" or "mostdense"
#'
#' @param k_value                 If method="mostdense," then a k value is required. Default is 400.
#'                                See get_cluster_max_densities for more details.
#'
#' @return A dataframe of solutions, one per cluster, in long format (one row per signature per cluster)
#'
#' @export

get_cluster_summary_points <- function(exposures_mcmc_output, exposure_mcmc_clusters, method = 'bestfit', k_value = 400) {
  if (method == 'mostdense') {
    max_density_per_cluster <- get_cluster_max_densities(
      exposures_mcmc_output, 
      exposure_mcmc_clusters,
      k_value
    )
  } else if (method == 'bestfit') {
    best_fit_per_cluster <- get_cluster_best_fits(exposures_mcmc_output, exposure_mcmc_clusters)
  } else {
    stop('Invalid value provided to argument "method"')
  }
}


