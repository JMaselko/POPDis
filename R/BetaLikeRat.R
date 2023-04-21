#' @title Beta Loglikelihood ratio
#'
#' @description This function assesses the goodness of fit of a given vector for two candidate Beta distributions. Calculates Beta likelihood ratio for two Beta distributions with two sets of parameters and a data vector.
#'
#' @param X data vector.
#' @param a1 shape parameter 1 of first candidate Beta distribution.
#' @param b1 shape parameter 2 of first candidate Beta distribution.
#' @param a2 shape parameter 1 of second candidate Beta distribution.
#' @param b2 shape parameter 2 of second candidate Beta distribution.
#' @param minDet value to replace any 0s or (1-minDet) for 1s, default 1e-9.
#' @return A statistic with large values indicating more support for the first candidate Beta distribution.
#' @examples
#' X=c(0.2, 0.4, 0.3, 0.6, 0.2, 0.1, 0.4)
#' BetaLikeRat(X, a1=1, b1=1, a2=10, b2=5)
#' @export

BetaLikeRat<-function(X,a1,b1,a2,b2,minDet=1e-9){
  n<-length(X)
  #Need to replace 0s and 1s:
  X<-replace(X,(X==0),minDet)
  X<-replace(X,(X==1),(1-minDet))
  bet1<-beta(a1,b1)
  if(bet1==0){bet1=minDet}
  bet2<-beta(a2,b2)
  if(bet2==0){bet2=minDet}
  #  bet<-n*(log(beta(a2,b2))-log(beta(a1,b1)))
  bet<-n*(log(bet2)-log(bet1))
  su<-sum((a1-a2)*log(X)+(b1-b2)*log(1-X))
  return(bet+su)
}
