library(nimble)
library(coda)
source("sim.RT.MVhyper.R")
source("NimbleModelRT MVhyper hurdleZTPois.R")
source("NimbleFunctionsRT MVhyper hurdleZTPois.R")
source("init.RT.MVhyper.R")
source("sSampler.R")

#If using Nimble version 0.13.1 and you must run this line 
nimbleOptions(determinePredictiveNodesInModel = FALSE)
# #If using Nimble before version 0.13.1, run this line instead
# nimble:::setNimbleOption('MCMCjointlySamplePredictiveBranches', FALSE) 

####Simulate some data####
N=50
#detection parameters
p0=0.25 #baseline detection probability
lambda=1 #count parameter given detection
sigma=0.5

#estimate of mean counts/capture event for lambda
mean(VGAM::rzapois(10000,lambda,pobs0=0))
#pmf for counts of 1:10
round(VGAM::dzapois(1:10,lambda,pobs0=0),2)

K=5 #number of occasions
buff=3 #state space buffer
X<- expand.grid(3:11,3:11) #make a trapping array
n.select=1 #samples to choose per trap/occasion

data=sim.RT.MVhyper(N=N,p0=p0,lambda=lambda,sigma=sigma,K=K,X=X,buff=buff,n.select=n.select,obstype="hurdleZTPois")

#What is the observed data?
str(data$y.ID) #the observed ID detections
head(data$this.j) #the trap of detection for all unidentified detections
head(data$this.k) #occasion of capture

#proportion samples ID'd
sum(data$y.ID)/sum(data$y.true)

##Fit model in Nimble##
#data augmentation level
M=150

J=nrow(X) #number of detectors
K2D=data$K2D #pull out trap operation

inits=list(p0=0.5,sigma=1,lambda=1) #ballpark inits to build data

#This function structures the simulated data to fit the model in Nimble (some more restructing below)
#Also checks some inits
nimbuild=init.RT.MVhyper(data,inits,M=M,obstype="hurdleZTPois")

#inits for nimble
Niminits <- list(z=nimbuild$z,s=nimbuild$s,ID=nimbuild$ID,capcounts=rowSums(nimbuild$y.true),
                 y.true=nimbuild$y.true,p0=inits$p0,sigma=inits$sigma,lambda=inits$lambda)

#constants for Nimble
J=nrow(data$X)
constants<-list(M=M,J=J,K=K,K2D=K2D,n.samples=nimbuild$n.samples,xlim=data$xlim,ylim=data$ylim)

# Supply data to Nimble. Note, y.true is completely latent.
z.data=c(rep(1,data$n.ID),rep(NA,M-data$n.ID))

Nimdata<-list(y.true=array(NA,dim=c(M,J,K)),
              ID=rep(NA,nimbuild$n.samples),z=z.data,X=as.matrix(X),capcounts=rep(NA,M))

# set parameters to monitor
parameters<-c('psi','p0','sigma','lambda','N','n')
nt=1 #thinning rate
#can also monitor a different set of parameters with a different thinning rate
parameters2 <- c("ID")
nt2=50#thin more

# Build the model, configure the mcmc, and compile
start.time<-Sys.time()
Rmodel <- nimbleModel(code=NimModel, constants=constants, data=Nimdata,check=FALSE,
                      inits=Niminits)
#can use "nodes" argument in configureMCMC below to omit y.true that is replaced below for faster
#configuration
conf <- configureMCMC(Rmodel,monitors=parameters, thin=nt,
                      monitors2=parameters2,thin2=nt2,useConjugacy = TRUE) 

#conf$printSamplers() #shows the samplers used for each parameter and latent variable

###One *required* sampler replacements

##Here, we remove the default samplers for y.true and y.event, which are not correct
#and replace it with the custom "IDSampler"
conf$removeSampler("y.true")
#precalculate these
n.jk=apply(nimbuild$y.true,c(2,3),sum) #
n.ID.jk=apply(data$y.ID,c(2,3),sum)
denom.choose=matrix(NA,J,K) #MVHyper denominator for all j,k
for(j in 1:J){
  for(k in 1:K){
    denom.choose[j,k]=NimChoose(n.jk[j,k],n.ID.jk[j,k])
  }
}

conf$addSampler(target = paste0("y.true[1:",M,",1:",J,",1:",K,"]"),
                type = 'IDSampler',control = list(M=M,J=J,K=K,K2D=K2D,n.samples=nimbuild$n.samples,
                                                  this.j=nimbuild$this.j,this.k=nimbuild$this.k,
                                                  y.ID=nimbuild$y.ID,n.jk=n.jk,denom.choose=denom.choose),
                silent = TRUE)

###Two *optional* sampler replacements:

#replace default activity center sampler that updates x and y locations separately with a joint update
#a little more efficient. sSampler below only tunes s when z=1. Should better tune activity centers for 
#uncaptured individuals. If this doesn't work well, try setting a fixed scale (adaptive=FALSE) in the range of 1/4 - 1/2 sigma.
conf$removeSampler(paste("s[1:",M,", 1:2]", sep=""))
for(i in 1:M){
  # conf$addSampler(target = paste("s[",i,", 1:2]", sep=""),
  #                 type = 'RW_block',control=list(adaptive=TRUE),silent = TRUE)
  conf$addSampler(target = paste("s[",i,", 1:2]", sep=""),
                  type = 'sSampler',control=list(i=i,xlim=nimbuild$xlim,ylim=nimbuild$ylim,
                                                 scale=1,adaptive=TRUE),silent = TRUE)
  #scale parameter here is just the starting scale. It will be tuned, unless adaptive=FALSE.
}

#replace independent lam0 and sigma samplers with block sampler better accommodating for posterior covariance
#should improve mixing and increase posterior effective sample size. AF_slice works better than block RW.
#Need to not use this update or modify it when using lam0 or sigma covariates.
conf$removeSampler(c("p0","sigma"))
conf$addSampler(target = c("p0","sigma"),
                  type = 'AF_slice',
                  control = list(adaptive=TRUE),
                  silent = TRUE)

# Build and compile
Rmcmc <- buildMCMC(conf)
# runMCMC(Rmcmc,niter=1) #this will run in R, used for better debugging
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

# Run the model.
start.time2<-Sys.time()
Cmcmc$run(2500,reset=FALSE) #short run for demonstration. can keep running this line to get more samples
end.time<-Sys.time()
end.time-start.time  # total time for compilation, replacing samplers, and fitting
end.time-start.time2 # post-compilation run time


mvSamples = as.matrix(Cmcmc$mvSamples)
plot(mcmc(mvSamples[2:nrow(mvSamples),]))

data$n.cap #true number of captured individuals

#look at ID posteriors. Not removing any burnin here...
mvSamples2 = as.matrix(Cmcmc$mvSamples2)

check.sample=1
#posterior prob this sample belongs to each individual number
round(table(mvSamples2[,check.sample])/(nrow(mvSamples2)-1),2)
#truth (for simulated data sets)
data$ID[check.sample]
#These should match data$ID for all individuals with identified captures, 1, ... n.ID
#individuals with ID>n.ID will have different numbers than data$ID (but samples with same true ID should tend to have same posterior ID)
data$n.ID
