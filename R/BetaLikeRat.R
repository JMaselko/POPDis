# Calculates Beta likelihood ratio for two Beta distributions with two sets of parameters
# And a data vector
# X: data vector [0,1]
# a1, b1: alpha and beta shape parameters of first candidate Beta distribution
# a2, b2: alpha and beta shape parameters of the second candidate Beta distribution
BetaLikeRat<-function(X,a1,b1,a2,b2,minDet=0.0001){
  n<-length(X)
  #Need to replace 0s and 1s:
  X<-replace(X,(X==0),minDet)
  X<-replace(X,(X==1),(1-minDet))
  #  return((beta(a1,b1)/beta(a2,b2))^n)*prod((X^(a1-a2)*(1-X)^(b1-b2)))
  #  return(prod((X^(a1-a2)*(1-X)^(b1-b2))))
  bet1<-beta(a1,b1)
  if(bet1==0){bet1=minDet}
  bet2<-beta(a2,b2)
  if(bet2==0){bet2=minDet}
  #  bet<-n*(log(beta(a2,b2))-log(beta(a1,b1)))
  bet<-n*(log(bet2)-log(bet1))
  su<-sum((a1-a2)*log(X)+(b1-b2)*log(1-X))
  return(bet+su)
}
