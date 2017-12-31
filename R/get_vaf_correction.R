#' Computes the VAF correction factor for a table of mutations
#'
#' The VAF correction factor adjusts the expected mean variant allele fraction
#' 
#' The correction factor is computed as 
#'
#' \deqn{\frac{p}{C^{(t)} p + C^{(n)} (1-p)}},
#'
#' where \emph{p} is the tumour content / purity, $C^{(t)}$ is the tumour
#' copy number, and $C^{(n)}$ is the normal copy number.
#'
#' @param mutation_table        The same mutation table used as input to
#'                              \code{\link{get_population_signatures}}.
#'
#' @return Numeric vector of correction factors

get_vaf_correction <- function(mutation_table) {
    with(
        mutation_table,
        tumour_content / ((tumour_content * tumour_copy) + ((1 - tumour_content) * normal_copy)),
    )
}
