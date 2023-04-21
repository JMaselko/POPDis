#' @title Time to Stationarity Slope
#'
#' @description Function calculates the significance of the slope for the mean pairwise Fst for a given dispersal model and time frame.
#'
#' @param Dispersal Use the specified Dispersion matrix, or alternatively if "Random" calls Dispers.Rand function to create random dispersal matrix at each year step.
#' @param Dispersion Specify the dispersal matrix
#' @param N Effective population size
#' @param PopLoc Number of populations
#' @param MinAge Minimum age for spawning cohort
#' @param MaxAge Maximum lifespan
#' @param loci Number of loci to simulate
#' @param Mutation Whether to incorporate mutation
#' @param Selection Whether to incorporate selection
#' @param selstr Selection strength
#' @param Stoch Whether to include stochasticity
#' @param mut Mutation rate
#' @param YStart Start year for recording stationarity variables
#' @param YInterval Interval years between testing slope
#' @param YFinal Stop simulation after these many generations
#' @param tolerance The number of significant digits to use for calculating significance of slope test

#' @return A data frame containing two year values
#' @returns SlopeMin year when first occurrence of significant 0 slope for mean and variance of Fst.
#' @returns SlopeMax year when last occurrence of significant 0 slope for mean and variance of Fst.
#'
#' @examples
#' Disp.StationSlope(Dispersal="Random", Dispersion=NA, N=1e3, PopLoc=10, MinAge=5, MaxAge=20, loci=1000, Mutation=TRUE, Selection=FALSE, selstr=0, Stoch=TRUE, mut=1e-6, YStart=1, YInterval=100, YFinal=1e3, tolerance=6)
#' @export
#' @import stats
#' @import smatr

Disp.StationSlope<-function(Dispersal="new", Dispersion=NA, N=1e3, PopLoc=10, MinAge=5, MaxAge=20, loci=1000, Mutation=TRUE, Selection=FALSE, selstr=0, Stoch=TRUE, mut=1e-6, YStart=1, YInterval=100, YFinal=1e3, tolerance=6){
  if (Dispersal=="Random"){ Dispersion=Dispers.Rand(Pops=PopLoc)}

  Disp.Stat<-data.frame(yr1=NA, meanFstSlope=NA, meanVarSlope=NA)#,
  VarMean<-data.frame(y=c(1:YInterval),meanFst=NA,varFst=NA )
  varcount<-1
  if (any(is.na(Dispersion))==FALSE) {PopLoc<-length(Dispersion[,1])}
  Pop=array(1,c(PopLoc,MaxAge,loci),dimnames=list(Population=c(1:PopLoc),Age=c(1:MaxAge),Loci=c(1:loci)))

  # Here I introduce the initial population allele frequencies as drawn from beta
  for (i in 1:PopLoc){
    Pop[i,1,]=rbeta(loci,0.5,0.5) # beta(0.5,0.5) is the non-informative Jeffrey's prior
    for (k in 1:MaxAge){Pop[i,k,]=Pop[i,1,] }
   }

  Ndisp<-matrix(data=0, nrow=PopLoc, ncol=PopLoc)
for (yr in 1:YFinal){
    if (Dispersal=="Random"){ Dispersion=Dispers.Rand(Pops=PopLoc)} #else

    ###########################
    #######################
    # Calculating Pop within function:
    # Step 1:
    # create allele frequencies for each spawning group at each location:
    x<-matrix(data=NA, nrow=PopLoc, ncol=loci)
    if (Selection==TRUE) {
      s<-runif(loci,-selstr,selstr)}

    for (k in 1:PopLoc){
      # Create spawning aggregates from each spawning cohort
      a_r<-0
      lR<-length(MinAge:MaxAge)#How many spawning cohorts
      for (r in MinAge:MaxAge){
        # Draw a sample of alelles from each spawner age group and combine:
        if (anyNA(Pop[k,r,])==TRUE){
          lR=lR-1 # Note this is a failure in cohort recruitment to spawning population
        } else {
          if (Stoch==TRUE){
            a_r <- a_r+rbinom(loci,2*N,prob=Pop[k,r,])  # stochastic number of alleles from each cohort
          } else { a_r <- a_r+round(2*N*Pop[k,r,]) }  # deterministic number of alleles from each cohort
        }
      }
      a_r<-a_r/(2*N*lR) # These are the allele frequencies in cohort aggregated spawner group
      # Now the Spawners of each population create babies:
      if (Mutation == TRUE) {
        a_r<- a_r*(1-mut)+(1-a_r)*mut # with mutation (p(major allele not mutating)+p(minor allele mutating))
      }     # No mutation
      # These are the larval allele frequencies from each source spawning location
      # where rows are populations and columns are loci:
      if (Selection==TRUE) {
        a_r<-((1+s)*a_r)/((1+s)*a_r+(1-s)*(1-a_r))
      }
      x[k,]<-a_r
      # creating dispersal matrix
      if (Stoch==TRUE){
        if (sum(Dispersion[k,]==0)){
          Ndisp[k,]=0
        } else {
          Ndisp[k,]<-rmultinom(1,N,Dispersion[k,])
        }
      }
    }
    # Below I create age 0 cohort at each of the destination locations
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
    #######################
    Fstmat=Fst.mat(Pop)
    VarMean[varcount,]<-c(yr,round(mean(Fstmat), tolerance),round(var(Fstmat), tolerance))

    if(varcount==YInterval){
      Disp.Stat<-rbind(Disp.Stat, c(yr, smatr::slope.test(VarMean$meanFst, VarMean$y)$p, smatr::slope.test(VarMean$varFst, VarMean$y)$p))
      }
    varcount<-varcount+1
    if(varcount>YInterval){varcount=1} # resetting so that I only look at the last YInterval number
  }

  Disp.Stat=na.omit(Disp.Stat)
  if (is.err(min(Disp.Stat$yr1[(which(Disp.Stat$meanFstSlope==1 & Disp.Stat$meanVarSlope==1))]))==TRUE){
    minSlopeYear<-YFinal
  } else {minSlopeYear<-min(Disp.Stat$yr1[(which(Disp.Stat$meanFstSlope==1 & Disp.Stat$meanVarSlope==1))])}

  maxSlopeYear<-max(na.omit(Disp.Stat$yr1[(which(Disp.Stat$meanFstSlope!=1 | Disp.Stat$meanVarSlope!=1))+1]))

  Result<-data.frame(SlopeMin=minSlopeYear, SlopeMax=maxSlopeYear)

  return(Result)
}
