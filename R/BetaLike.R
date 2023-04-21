#' @title Beta Loglikelihood
#'
#' @description This function calculates the log likelihood of a given vector given the alpha and beta shape parameters of the Beta distribution. It uses the Stirling's approximation of the beta function.
#`
#' @param X data vector.
#' @param a1 shape parameter 1 of the Beta distribution.
#' @param b1 shape parameter 2 of the Beta distribution.
#' @param minDet value to replace any 0s or (1-minDet) for 1s, default minDet=1e-9.
#' @return A statistic.
#' @examples
#' X=c(0.2, 0.4, 0.3, 0.6, 0.2, 0.1, 0.4)
#' BetaLike(X, a1=1, b1=1)
#'
#' @export
#' @import stats

BetaLike<-function(X,a1,b1,minDet=1e-9){
  n<-length(X)
  pi<-base::pi
  #Need to replace 0s and 1s:
  X<-replace(X,(X==0),minDet)
  X<-replace(X,(X==1),(1-minDet))
  # Using the Stirling's approximation of beta(a)
  b<-log(sqrt(2*pi))+(a1-1/2)*log(a1)+(b1-1/2)*log(b1)-(a1+b1-1/2)*log(a1+b1)
  bet<-(-n*b)
  su<-sum((a1-1)*log(X)+(b1-1)*log(1-X))
  return(bet+su)
}
