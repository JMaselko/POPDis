# Beta Loglikelihood.
# Calculates the log likelihood of a given vector
# and alpha and beta parameters of the beta distribution.
# It uses the Stirling's approximation of the beta function

BetaLike<-function(X,a1,b1,minDet=0.0001){
  n<-length(X)
  pi<-base:::pi
  #Need to replace 0s and 1s:
  X<-replace(X,(X==0),minDet)
  X<-replace(X,(X==1),(1-minDet))
  # Using the Stirling's approximation of beta(a)
  b=log(sqrt(2*pi))+(a1-1/2)*log(a1)+(b1-1/2)*log(b1)-(a1+b1-1/2)*log(a1+b1)
  bet=-n*b
  su<-sum((a1-1)*log(X)+(b1-1)*log(1-X))
  return(bet+su)
}
