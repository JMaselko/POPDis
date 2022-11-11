# Calculates number of time steps until stationarity is reached.
# Inputs are:
# Disp.matrix: Dispersal matrix
# delta: tolerance, i.e. stop when the change in matrices is less than the tolerance parameter
# steps: maximum number of iterations
# This is based on x%*%x
########################
TimetoStationarity<-function(Disp.matrix=c(), delta=1e-6,steps=1e6){

  x<-Disp.matrix
  k=nrow(x)
  x<-x/rowSums(x)[row(x)]
  x[is.na(x)]<-0
  x
  #delta = 1e-6 # tolerance to denote that we reached stationariry

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
