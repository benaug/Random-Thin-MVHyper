NimModel <- nimbleCode({
  #detection function priors
  p0 ~ dunif(0,1)
  sigma ~ dunif(0,20)
  #sample deposition rate prior
  lambda ~ dunif(0,20)
  #data augmentation prior
  psi ~ dunif(0,1)
  #likelihoods (except for s priors)
  for(i in 1:M) {
    z[i] ~ dbern(psi)
    s[i,1] ~ dunif(xlim[1],xlim[2])
    s[i,2] ~ dunif(ylim[1],ylim[2])
    pd[i,1:J] <- GetDetectionProb(s = s[i,1:2], X = X[1:J,1:2], J=J,sigma=sigma, p0=p0, z=z[i])
    y.true[i,1:J,1:K] ~ dHurdleZTPoisVector(pd=pd[i,1:J],K2D=K2D[1:J,1:K],z=z[i],lambda=lambda) #vectorized obs mod
  }
  #derived variables
  capcounts[1:M] <- Getcapcounts(ID=ID[1:n.samples],M=M,capcounts.ID=capcounts.ID[1:M]) #intermediate object
  n <- Getncap(capcounts=capcounts[1:M]) #number of captured individuals
  N <- sum(z[1:M])
})# end model
