# Function to test for number of generations to stationarity by calculating correlation, mean and variance for allele frequencies between YInterval generations:
# Input parameters are:
# Dispersal="new": Use "Random" or "Panmix" to generate corresponding dispersal matrices
# Dispersion=NA: Specify the dispersal matrix
# N=1e3: Effective population size
# PopLoc=10: Number of populations
# MinAge=5: Minimum age for spawning cohort
# MaxAge=20: Maximum lifespan
# loci=1000: Number of loci to simulate
# Mutation=TRUE: Wheteher to incorporate mutation
# Selection=FALSE: Whether to incorporate selection
# selstr=0: Selection strength
# Stoch=TRUE: Whether to include stochasticity
# mut=1e-6: Mutation rate
# YStart=1: Start year for recording stationarity variables
# YInterval=10: Interval years between recording stationarity variables
# YFinal=1e3: Stop simulation after these many generations

Disp.Stationarity.sim<-function(Dispersal="new",Dispersion=NA, N=1e3,PopLoc=10,MinAge=5, MaxAge=20,loci=1000, Mutation=TRUE,Selection=FALSE, selstr=1,Stoch=TRUE, mut=1e-6,YStart=1,YInterval=10,YFinal=1e3){
  #   set.seed(1)
  if (Dispersal=="Random"){ Dispersion=Dispers.Rand(Pops=PopLoc)} #else
  if (Dispersal=="Panmix"){ Dispersion=Dispers.Panmix(Pops=PopLoc)}
  Disp.Stat<-data.frame(yr1=NA,yr2=NA,corr=NA,meanFst=NA, varFst=NA)
  if (any(is.na(Dispersion))==FALSE) {PopLoc<-length(Dispersion[,1])}
  Pop=array(1,c(PopLoc,MaxAge,loci),dimnames=list(Population=c(1:PopLoc),Age=c(1:MaxAge),Loci=c(1:loci)))

  # Here I introduce the initial population allele frequencies as drawn from beta
  for (i in 1:PopLoc){
    Pop[i,1,]=rbeta(loci,0.5,0.5) # beta(0.5,0.5) is the non-informative Jeffrey's prior
    for (k in 1:MaxAge){Pop[i,k,]=Pop[i,1,] }
  }
  i=1
  yr1<-1
  F.mat<-c()
  F.mat1<-c()
  Ndisp<-matrix(data=0, nrow=PopLoc, ncol=PopLoc)
  for (yr in 1:YFinal){
    # calling user defined function:
    if (Dispersal=="Random"){ Dispersion=Dispers.Rand(Pops=PopLoc)}
    # Calculating Pop within function:
    # Step 1:
    # create allele frequencies for each spawning group at each location:
    x<-matrix(data=NA, nrow=PopLoc, ncol=loci)
    if (Selection==TRUE) {
      s<-runif(loci,-selstr,selstr)}
    Nc<-N
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
            a_r <- a_r+rbinom(loci,2*Nc,prob=Pop[k,r,])  # stochastic number of alleles from each cohort
          } else { a_r <- a_r+round(2*Nc*Pop[k,r,]) }  # deterministic number of alleles from each cohort
        }
      }
      a_r<-a_r/(2*Nc*lR) # These are the allele frequencies in cohort aggregated spawner group
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
        if (sum(Dispersion[k,])==0){
          Ndisp[k,]=0
        } else {
          Ndisp[k,]<-rmultinom(1,N,Dispersion[k,])
        }
       }
    }
    # Below I create age 0 cohort at each of the destination locations
    x_0<-matrix(data=0, nrow=PopLoc, ncol=loci)
    ### Now for a more streamlined approach:
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

    if(yr %in% seq(YStart,YFinal, by=YInterval)){
      F.mat<-Fst.mat(Pop)
      if(yr>YStart){
        if (is.err(cor(F.mat,F.mat1))){corF=0} else {corF=cor(F.mat,F.mat1)}
        Disp.Stat<-rbind(Disp.Stat, c(yr1,yr,corF, mean(F.mat), var(F.mat)))
      }
      yr1<-yr
      F.mat1<-F.mat
      i=i+1
    }
  }
  return(Disp.Stat)
}
