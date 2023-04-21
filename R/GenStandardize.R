#' @title Genotype Standardization
#'
#' @description This function standardizes genotypes for PCA based on Coop et al. 2010.
#`
#' @param Gen Genotype matrix of n individual by l loci
#' @return A matrix.
#' @examples
#' Gen=matrix(rbinom(100, 2, 0.5), ncol=20)
#' GenStandardize(Gen)
#' @export

GenStandardize<-function(Gen){
  if (any(apply(Gen,2,var)==0)){  Gen=Gen[,-which(apply(Gen,2,var)==0)]}
  pl<-colSums(Gen)/(2*nrow(Gen))
  den<-1/sqrt(2*pl*(1-pl))
  num<-t(apply(Gen,1,'-',2*pl))
  res=num*den
  return(res)
}
