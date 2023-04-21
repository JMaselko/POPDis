#' @title Dispersal Generated Allele Frequencies
#'
#' @description This function the main stochastic genetic population algorithm that generates the allele frequencies based on dispersal matrix.
#'
#' @param Dispersal  Use the specified Dispersion matrix, or alternatively if "Random" calls Dispers.Rand function to create random dispersal matrix at each year step.
#' @param Dispersion  Matrix of dispersion probabilities.
#' @param N Effective population size
#' @param PopLoc Number of populations to model
#' @param YearSteps  Number of year steps to simulate
#' @param MinAge Minimum age of breeding cohort
#' @param MaxAge Maximum age of oldest cohort, i.e. maximum lifespan
#' @param loci Number of loci to simulate
#' @param Mutation Whether to include a mutation parameter in the simulation
#' @param Selection Whether to include a selection parameter in the simulation
#' @param selstr Selection strength parameter [-1,1] where each cohort selection~unif(1, -selstr, selstr)
#' @param selrand Whether selection is random for each loci~unif(loci, -selstr, selstr)
#' @param Stoch Whether to model the dispersion stochastically where each generation parents are chosen at random
#' @param mut mutation rate
#' @return A list containing allele frequencies for each age class and population
#' @examples
#' Disp.sim(Dispersal="Random", Dispersion=NA, N=1e3, PopLoc=10, YearSteps=1000, MinAge=5, MaxAge=20,loci=1000, Mutation=TRUE, Selection=FALSE, selstr=0, selrand=FALSE, SStoch=TRUE, mut=1e-6)
#' @export
#'

Disp.sim<-function(Dispersal="new", Dispersion=NA, N=1e3, PopLoc=10, YearSteps=1000, MinAge=5, MaxAge=20, loci=1000, Mutation=TRUE, Selection=FALSE, selstr=0, selrand=FALSE, Stoch=TRUE, mut=1e-6){

   # If Dispersion Matrix is not provided:
  if (Dispersal=="Random"){ Dispersion=Dispers.Rand(Pops=PopLoc)}
  if (any(is.na(Dispersion))==FALSE) {PopLoc<-length(Dispersion[,1])}

  Pop=array(1,c(PopLoc,MaxAge,loci),dimnames=list(Population=c(1:PopLoc),Age=c(1:MaxAge),Loci=c(1:loci)))

  # Here I introduce the initial population allele frequencies as drawn from beta
  for (i in 1:PopLoc){
    Pop[i,1,]=rbeta(loci,0.5,0.5) # beta(0.5,0.5) is the non-informative Jeffrey's prior
    for (k in 1:MaxAge){Pop[i,k,]=Pop[i,1,] }
   }
  for (t in 1:YearSteps){
    # Randomly choose dispersion parameters simulating panmixia:
    if (Dispersal=="Random"){Dispersion=Dispers.Rand(Pops=PopLoc)}
    # Step 1:
    # create allele frequencies for each spawning group at each location:
    x<-matrix(data=NA, nrow=PopLoc, ncol=loci)
    Ndisp<-matrix(data=0, nrow=PopLoc, ncol=PopLoc)
    # Cohort specific selection
    if (Selection==TRUE){
      if (selrand==TRUE){
        s<-runif(loci,-selstr,selstr) # when loci experience independent selection (up to selstr) and direction
        } else if (selrand==FALSE){
        s<-rep(runif(1,-selstr,selstr), loci)} # when all loci experience equal selection strength (up to selstr) and direction
    }
    for (k in 1:PopLoc){
      # Create spawning aggregates from each spawning cohort
      a_r<-0
      lR<-length(MinAge:MaxAge)# Number of spawning cohorts
      for (r in MinAge:MaxAge){
        # Draw a sample of alelles from each spawner age group and combine:
        if (anyNA(Pop[k,r,])==TRUE){ # Note this is a failure in cohort recruitment to spawning population
          lR=lR-1
        } else {
          if (Stoch==TRUE){
            a_r <- a_r+rbinom(loci,2*N,prob=Pop[k,r,])  # stochastic number of alleles from each cohort
          } else { a_r <- a_r+round(2*N*Pop[k,r,]) }  # deterministic number of alleles from each cohort
        }
      }
      a_r<-a_r/(2*N*lR) # These are the allele frequencies in cohort aggregated spawner group
      # Now the aggregate cohort Spawners of each population create babies:
      if (Mutation == TRUE) {
        a_r<- a_r*(1-mut)+(1-a_r)*mut # with mutation (p(major allele not mutating)+p(minor allele mutating))
      }     # No mutation
      # The a_r are the allele frequencies in the larvae that are then subject to selection
      # during their simulated pelagic larval stage:
      if (Selection==TRUE) {
        a_r<-((1+s)*a_r)/((1+s)*a_r+(1-s)*(1-a_r))
      }
      # These are the larval allele frequencies from each source spawning location
      # where rows are populations and columns are loci:
      x[k,]<-a_r
      # creating dispersal matrix
      if (Stoch==TRUE){
        if (sum(Dispersion[k,])==0){
          Ndisp[k,]=0
        } else {
          Ndisp[k,]<-rmultinom(1,N,Dispersion[k,])
        }
      }
    }
    # Below I create age 0 cohort that is joining the destination population based
    # on the Dispersion matrix:
    x_0<-matrix(data=0, nrow=PopLoc, ncol=loci)

    if (Stoch==FALSE){
      Ndisp=as.matrix(round(Dispersion*N))
    }
    x_0<-t(Ndisp/colSums(Ndisp)[col(Ndisp)])%*%x
    if (anyNA(colSums(x_0))){
      x_0[,which(is.na(colSums(x_0)))]<-x[,which(is.na(colSums(x_0)))] # this is when failure of any migrants
    }
     # Drop last cohort, and update each cohort by 1 and add year 1 from x_0
    if (MaxAge==1){
      Pop[,1,]<-x_0
    } else {
      for (c in MaxAge:2){
        Pop[,c,]=Pop[,(c-1),]
      }
      Pop[,1,]<-x_0
    }
  }
  return(Pop)
}
