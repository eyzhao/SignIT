% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_mutation_types.R
\name{collapse_mutation_types}
\alias{collapse_mutation_types}
\title{Collapses mutation types of G>N and A>N mutations into C>N and T>N mutations
so that there are 6 possible base changes.}
\usage{
collapse_mutation_types(mutation_types)
}
\arguments{
\item{mutation_types}{Character vector of mutation types in format C[T>A]G
as an example of a T>A mutation in CTG trinucleotide context.}
}
\value{
New mutation types with any G>N/A>N mutations replaced and their
        trinucleotide contexts reverse complemented
}
\description{
Collapses mutation types of G>N and A>N mutations into C>N and T>N mutations
so that there are 6 possible base changes.
}
