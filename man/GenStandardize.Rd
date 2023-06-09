% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GenStandardize.R
\name{GenStandardize}
\alias{GenStandardize}
\title{Genotype Standardization}
\usage{
GenStandardize(Gen)
}
\arguments{
\item{Gen}{Genotype matrix of n individual by l loci}
}
\value{
A matrix.
}
\description{
This function standardizes genotypes for PCA based on Coop et al. 2010.
}
\examples{
Gen=matrix(rbinom(100, 2, 0.5), ncol=20)
GenStandardize(Gen)
}
