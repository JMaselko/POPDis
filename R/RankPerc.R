#' @title Rank Percentiles
#'
#' @description This function rescales the dispersal values for for graphing scaling visualization of relative proportions.
#`
#' @param x matrix of values
#' @return A matrix.
#' @examples
#' x=matrix(runif(100,min=0, max=1), ncol=20)
#' RankPerc(x)
#' @export

RankPerc<-function(x){
  y<-rank(x)
  y=(y-min(y))/(max(y)-min(y))
  return(matrix(ncol=ncol(x),nrow=nrow(x), byrow=F, data=y))
}
