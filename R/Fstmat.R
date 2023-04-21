#' @title Pairwise Gst matrix
#'
#' @description This function calculates the pairwise G"st matrix between K populations averaging over age classes (see Nei and Chesser, 1983)
#`
#' @param Pop Population allele frequency object as generated from Disp.sim function
#' @return A distance matrix.
#' @examples
#' Pop=Disp.sim(Dispersal="Random", Dispersion=NA, N=1e3, PopLoc=10, YearSteps=1000, MinAge=5, MaxAge=20,loci=1000, Mutation=TRUE, Selection=FALSE, selstr=1, Stoch=TRUE, mut=1e-6)
#' Fst.mat(Pop)
#'
#' @export

Fst.mat<-function(Pop){
  PopLoc<-dim(Pop)[1]
  FstMat=matrix(nrow=PopLoc,ncol=PopLoc)
  for (i in 1:PopLoc){
    for (j in 1:PopLoc){
      if(length(Pop[1,,1])==1){
        FstMat[i,j]=stGst(Pop1=(na.omit(Pop[i,,])),Pop2=(na.omit(Pop[j,,])))
      } else {
        FstMat[i,j]=stGst(Pop1=colMeans(na.omit(Pop[i,,])),Pop2=colMeans(na.omit(Pop[j,,])))
      }
    }
  }

 return(as.dist(FstMat))
}
