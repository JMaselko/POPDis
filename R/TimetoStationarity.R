#' @title Matrix Stationarity
#'
#' @description This function calculates time to stationarity based on the x%*%x' approach
#`
#' @param Disp.matrix Dispersal matrix
#' @param delta Tolerance parameter to denote stationarity
#' @param steps maximum number of iterations
#' @return A distance matrix.
#' @examples
#' Disp.matrix=matrix(data=runif(100,min=0,max=1), nrow=10)#' Fst.mat(Pop)
#' TimetoStationarity(Disp.matrix)
#' @export

TimetoStationarity<-function(Disp.matrix=c(), delta=1e-6,steps=1e6){
  x<-Disp.matrix
  k=nrow(x)
  x<-x/rowSums(x)[row(x)]
  x[is.na(x)]<-0
  x
  delta.obs = delta+1 # just something larger to start the while loop

  num.steps = steps # number of steps to run the markov chain

  current.matrix = vector(mode= "list") # store the matrices through time steps of simulation

  current.matrix[[1]] = x

  for(i in 2:num.steps){
    current.matrix[[i]] = current.matrix[[i-1]]%*%x
    delta.obs = sum(abs(current.matrix[[i]]-current.matrix[[i-1]]))
    if(delta.obs<delta){

      break
    }
    if(delta.obs>delta & i == num.steps){

    }
  }
  return(i)
}
