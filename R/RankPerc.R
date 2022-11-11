# Rank Percentage for better heatmap
# Converts the vector values to percentiles of ranks
# It results in better graphing for heatmaps
RankPerc<-function(x){
  y<-rank(x)
  y<-normalize(y/max(y))
  return(matrix(ncol=ncol(x),nrow=nrow(x), byrow=F, data=y))
}
