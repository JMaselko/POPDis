# Estimate Beta distribution shape parameters alpha and beta
# from mean and variance using method of moments
BetaParEst <- function(m, v) {
  a <- (1/v)*((1 - m)*m^2) - m
  #  b <- (m/v)*(1-m)^2-1+m
  b <- a/m-a
  return(params = list(alpha = a, beta = b))
}
