#' @title Random Dispersal
#'
#' @description This function calculates the rescaled dispersal matrix with rows summing to 1
#'
#' @param Pops Number of populations. Default Pops=10.
#' @return A Pops*Pops matrix of dispersal parameters.
#' @examples
#' Dispers.Rand(Pops=10)
#' @export

Dispers.Rand<-function(Pops=10){
  Dispers.Rand=matrix(runif(Pops*Pops,max=1/(Pops)),nrow=Pops,ncol=Pops)
  Dispers.Rand=Dispers.Rand/rowSums(Dispers.Rand)[row(Dispers.Rand)]
  return(Dispers.Rand)
}
