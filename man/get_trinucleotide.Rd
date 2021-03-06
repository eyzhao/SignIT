% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_mutation_types.R
\name{get_trinucleotide}
\alias{get_trinucleotide}
\title{Gets the trinucleotide context for a set of genomic coordinates.}
\usage{
get_trinucleotide(chr, pos, ref, genome = NULL)
}
\arguments{
\item{chr}{Character vector of chromosome names. Must match provided genome.}

\item{pos}{Integer vector of genomic positions.}

\item{genome}{BSgenome object - default is BSgenome.Hsapiens.UCSC.hg19}
}
\value{
Character vector of trinucleotide contexts
}
\description{
Gets the trinucleotide context for a set of genomic coordinates.
}
