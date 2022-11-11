# This function calculates the pairwise G"st value for two genotype vectors.
# It is adapted from:
# Meirmans and Hedricks 2011:Assessing population structure: FST and related measures
stGst=function(Pop1,Pop2){
  Hs<-1-mean(Pop1^2+Pop2^2)/2
  Ht<-mean(1-((Pop1+Pop2)/2)^2)
  #G"st formulation:
  sttGst<-2*(Ht-Hs)/((2*Ht-Hs)*(1-Hs))

  return(sttGst)
}
