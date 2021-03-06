% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_nnls_solution.R
\name{plot_nnls_solution}
\alias{plot_nnls_solution}
\title{Plots a simple NNLS solution}
\usage{
plot_nnls_solution(exposures_mcmc_output, signature_trim = "Signature")
}
\arguments{
\item{exposures_mcmc_output}{}
}
\value{
A ggplot object showing the NNLS solution for each exposure value.
}
\description{
Obtains the naive NNLS signature solution and plots the exposures.
}
\details{
Non-negative least squares computes the value of x which minimizes the residual sum
of squares \eqn{\lVert Ax = b \rVert^2} given a known matrix A and vector b.
This provides a simple solution to the mutation signature problem. However,
merely minimizing the residual sum of squares may result in overfitting the solution.
This plotting function computes and plots a simple NNLS-based solution to the
exposure vector and returns a plot. This can serve as a baseline comparison or sanity check
against the MCMC methods featured in SignIT.
}
