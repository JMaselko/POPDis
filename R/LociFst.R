#' @title Pairwise Fst for loci
#'
#' @description This function calculates the pairwsie Fst for each loci between two populations allele frequency matrices.
#`
#' @param Pop1 Population 1 allele frequency matrix
#' @param Pop2 Population 2 allele frequency matrix
#' @return A vector.
#' @examples
#' Pop1=matrix(rbeta(10,1,1))
#' Pop2=matrix(rbeta(10,1,10))
#' LociFst(Pop1, Pop2)
#'
#' @export

LociFst=function(Pop1,Pop2){
  He1<-(1-(Pop1^2+(1-Pop1)^2))
  He2<-(1-(Pop2^2+(1-Pop2)^2))
  Hs<-(He1+He2)/2
  pb<-(Pop1+Pop2)/2
  Ht<-2*pb*(1-pb)
  MFst<-((Ht)-(Hs))/(Ht)
  return(MFst)
}
