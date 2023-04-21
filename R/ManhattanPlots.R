#' @title Manhattan Plot
#'
#' @description This function creates a vector of pairwise Fst calculated for each loci between two populations i and j from POP object
#`
#' @param Pop Object created by Disp.sim function
#' @param i which 1st populations to compare
#' @param j which 2nd populations to compare
#' @param AllAg Whether to calculate allele frequencies across all age groups
#' @param Age which age to use for single age allele frequencies
#' @return A vector.
#' @examples
#' Pop=Disp.sim(Dispersal="Random", Dispersion=NA, N=1e3, PopLoc=10, YearSteps=1000, MinAge=1, MaxAge=1,loci=100, Mutation=TRUE, Selection=FALSE, selstr=1, Stoch=TRUE, mut=1e-6)
#' ManhattanPlots(Pop, i=1, j=2, AllAg=FALSE, Age=1)
#'
#' @export

ManhattanPlots<-function(Pop, i=1, j=2, AllAg=FALSE, Age=1){
 if (AllAg==TRUE){
  } else {
    Pop1=(na.omit(Pop[i,Age,]))
    Pop2=(na.omit(Pop[j,Age,]))
  }
  Fst<-LociFst(Pop1,Pop2)
  return(Fst)
}
