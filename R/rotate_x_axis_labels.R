#' Rotates the x axis labels 90 degrees in ggplot graphs.
#'
#' Convenience function for plots.
#'
#' @return ggplot component which can be added to a ggplot graph with the "+" operator.
#'
#' @import ggplot2

rotate_x_axis_labels <- function() {
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4))
}
