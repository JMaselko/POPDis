#' @title HWE test
#'
#' @description This function calculate HWE Chisquare test p-value for genotype matrix:
#`
#' @param Gen Genotype matrix of n individual by l loci
#' @return A vector.
#' @examples
#' Gen=matrix(rbinom(100, 2, 0.5), ncol=20)
#' HWE.chisq(Gen)
#' @export

HWE.chisq<-function(Gen){
  HWE.p<-c()
  n<-nrow(Gen)
  if (is.null(n)==FALSE){
    for (i in 1:ncol(Gen)){
      AA<-sum(Gen[,i]==0)
      Aa<-sum(Gen[,i]==1)
      aa<-sum(Gen[,i]==2)
      Obs.f<-(2*aa+Aa)/(2*n)
      p<-dbinom(0:2,2,Obs.f)
      chisq<-sum((c(AA,Aa,aa)-p*n)^2/(n*p))
      HWE.p[i]<-pchisq(chisq,1,lower.tail=FALSE)
    }
  }
  return(HWE.p)
}
