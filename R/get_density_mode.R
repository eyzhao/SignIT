#' Convenience function for estimating the mode of dataset's density using kernel estimation.
#' The kernel is estimated using R's base density() function and the max density point
#' is chosen. If multiple points have equally high estimated density, the central one is chosen.
#'
#' @param v         Vector of distributed data as input for kernel density estimation
#'
#' @return A value specifying the position of the density mode
#'
#' @import magrittr

get_density_mode <- function(v) {
    if (all(length(unique(v)) == 1)) {
        return(unique(v))
    } else {
        d <- density(v)
        return(d$x[d$y == max(d$y)] %>% .[ceiling(length(.) / 2)])
    }
}
