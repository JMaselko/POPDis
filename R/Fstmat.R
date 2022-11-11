# Calculating Pairwise Gst matrix
# It uses the stGst function:
# from the Pop list object:

Fst.mat<-function(Pop){
  PopLoc<-dim(Pop)[1]
  FstMat=matrix(nrow=PopLoc,ncol=PopLoc)
  for (i in 1:PopLoc){
    for (j in 1:PopLoc){
      if(length(Pop[1,,1])==1){
        FstMat[i,j]=stGst(Pop1=(na.omit(Pop[i,,])),Pop2=(na.omit(Pop[j,,]))) # MeanFst(Pop1=(na.omit(Pop[i,,])),Pop2=(na.omit(Pop[j,,])))
     } else {
        FstMat[i,j]=stGst(Pop1=colMeans(na.omit(Pop[i,,])),Pop2=colMeans(na.omit(Pop[j,,]))) # MeanFst(Pop1=colMeans(na.omit(Pop[i,,])),Pop2=colMeans(na.omit(Pop[j,,])))#
      }
    }
  }

 return(as.dist(FstMat))
}
