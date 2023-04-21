#' @title Vector Normalization
#'
#' @description This function normalizes a set of values to range (0:1). It is used in the RankPerc function for heatmap graphing.
#`
#' @param x vector of values
#' @return A vector.
#' @examples
#' x=runif(100,min=0, max=1)
#' normalize(x)
#' @export

normalize<-function(x){
  return((x-min(x))/(max(x)-min(x)))
}
