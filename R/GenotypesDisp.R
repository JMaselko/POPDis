# Creating n Genotypes from the Pop matrix:
# Inputs are:
# POP: Population allele frequency object created by Disp.sim function
# sampn: Vector of the number of individuals to sample from each population
# Allages=TRUE: whether to sample genotypes from all age groups or only a single age class
# Ag=1: Age cohort to sample

Genotypes.Disp<-function(POP=Disp.tst,sampn=rep(10,10), AllAges=TRUE, Ag=1){
  #sampn is number of individuals to sample per populations
  Populations<-length(POP[,1,1])
  loci<-length(POP[1,1,])
  ages<-length(POP[1,,1])
  Pop.genotypes<-list()
  r=1
  for (Pops in 1:Populations){
    if (sampn[Pops]>0){
      a=1
      b=1
      gen.matrix<-data.frame(matrix(NA,ncol=(1+loci),nrow=sampn[Pops]))#matrix(ncol=(1+loci),nrow=sampn[Pops])
      for (ns in 1:(sampn[Pops]/ages)){
        if (AllAges==FALSE){
          # sample genotypes from single age group Ag
          # genotypes are 0 (aa), 1 (Aa), 2 (AA)
          for (i in 1:ages){
            gen.matrix[b,]<-c(Ag,rbinom(loci,2,POP[Pops,Ag,]))
            b=b+1
          }
        } else {
          # sample genotypes randomly across all ages
          for (ag in 1:ages){
            gen.matrix[a,]<-c(ag,rbinom(loci,2,POP[Pops,ag,]))
            a<-a+1
          }
        }
      }
      colnames(gen.matrix)=c("Age",(1:loci))
      Pop.genotypes[[r]]<-as.matrix(gen.matrix)#na.omit(gen.matrix)
      r=r+1
    }
  }
  names(Pop.genotypes)=which(sampn>0)
  return(Pop.genotypes)
}
