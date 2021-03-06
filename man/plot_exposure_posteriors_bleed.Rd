% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_exposure_posteriors_bleed.R
\name{plot_exposure_posteriors_bleed}
\alias{plot_exposure_posteriors_bleed}
\title{Plot the signature exposure posterior distribution along with a visualization of
signature bleed interactions.}
\usage{
plot_exposure_posteriors_bleed(
  exposures_mcmc_output,
  view = "violin",
  min_bleed = 0.2,
  signature_trim = "Signature",
  legend = TRUE,
  plot_heights = c(1, 1),
  units = "mutations"
)
}
\arguments{
\item{exposures_mcmc_output}{Output from get_exposures.}

\item{view}{Can be either "violin" or "boxplot" (default: "violin")}

\item{min_bleed}{Bleed threshold, between 0 and 1. Default: 0.2.}

\item{signature_trim}{Text to be trimmed off of the x-axis signature
labels.}

\item{legend}{Boolean - whether or not to include the legend.}

\item{plot_heights}{Vector of 2 numbers denoting relative heights
of exposure plot and bleed graph.}

\item{units}{Units to present exposures in. Can be
'mutations', 'megabase', or 'fraction'.}
}
\value{
A cowplot merge of two ggplots showing the posterior distributions and signature bleed.
}
\description{
Plot the signature exposure posterior distribution along with a visualization of
signature bleed interactions.
}
