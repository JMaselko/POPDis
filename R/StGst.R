#' @title Meirmans and Hedricks Gst
#'
#' @description This function calculates the pairwise G"st value for two allele frequency vectors based on Meirmans and Hedricks 2011: "Assessing population structure: FST and related measures".
#`
#' @param Pop1 Population 1 allele frequency vector
#' @param Pop2 Population 2 allele frequency vector
#' @return A value.
#' @examples
#' Pop1=runif(100, min=0, max=1)
#' Pop2=runif(100, min=0, max=1)
#' stGst(Pop1, Pop2)
#' @export

stGst=function(Pop1,Pop2){
  Hs<-1-mean(Pop1^2+Pop2^2)/2
  Ht<-mean(1-((Pop1+Pop2)/2)^2)
  #G"st formulation:
  sttGst<-2*(Ht-Hs)/((2*Ht-Hs)*(1-Hs))

  return(sttGst)
}
