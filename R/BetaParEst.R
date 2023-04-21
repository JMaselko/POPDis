#' @title Beta distribution
#'
#' @description This function calculates the alpha and beta shape parameters of the Beta distribution via the method of moments.
#'
#' @param m data mean.
#' @param v data variance.
#' @return A list containing alpha and beta parameters of the Beta distribution.
#' @examples
#' BetaParEst(m=0.5, v=0.2)
#' @export

BetaParEst <- function(m, v) {
  a <- (1/v)*((1 - m)*m^2) - m
  #  b <- (m/v)*(1-m)^2-1+m
  b <- a/m-a
  return(params = list(alpha = a, beta = b))
}
