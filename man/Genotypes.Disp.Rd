% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GenotypesDisp.R
\name{Genotypes.Disp}
\alias{Genotypes.Disp}
\title{Genotype Generator}
\usage{
Genotypes.Disp(POP, sampn = rep(10, 10), AllAges = TRUE, Ag = 1)
}
\arguments{
\item{POP}{Population allele frequency object as generated from Disp.sim function}

\item{sampn}{Vector of the number of individuals to sample from each population}

\item{AllAges}{whether to sample genotypes from all age groups or only a single age class}

\item{Ag}{Age cohort to sample}
}
\value{
A list of genotypes with ages and populations.
}
\description{
This function generates the individual genotypes by drawing from Disp.sim allele frequency object
}
\examples{
Pop=Disp.sim(Dispersal="Random", Dispersion=NA, N=1e3, PopLoc=10, YearSteps=1000, MinAge=5, MaxAge=20,loci=1000, Mutation=TRUE, Selection=FALSE, selstr=1, Stoch=TRUE, mut=1e-6)
Genotypes.Disp<-function(Pop, sampn=rep(10,10), AllAges=TRUE, Ag=1)
}
