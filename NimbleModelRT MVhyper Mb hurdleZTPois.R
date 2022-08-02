NimModel <- nimbleCode({
  #detection function priors
  p0.p~dunif(0,1)
  p0.c~dunif(0,1)
  sigma~dunif(0,20)
  #sample deposition rate prior
  lambda~dunif(0,20)
  #data augmentation prior
  psi~dunif(0,1)
  #likelihoods (except for s priors)
  for(i in 1:M){
    z[i] ~ dbern(psi)
    s[i,1] ~ dunif(xlim[1],xlim[2])
    s[i,2] ~ dunif(ylim[1],ylim[2])
    kern[i,1:J] <- GetKern(s = s[i,1:2], X = X[1:J,1:2], J=J,sigma=sigma, z=z[i])
    pd.p[i,1:J] <- Getpd(kern = kern[i,1:J],p0=p0.p,J=J,z=z[i])
    pd.c[i,1:J] <- Getpd(kern = kern[i,1:J],p0=p0.c,J=J,z=z[i])
    y.true[i,1:J,1:K] ~ dHurdleZTPoisVector(pd.p=pd.p[i,1:J],pd.c=pd.c[i,1:J],y.state[i,1:J,1:K],
                                            K2D=K2D[1:J,1:K],z=z[i],lambda=lambda) #vectorized obs mod
  }
  #derived variables
  capcounts[1:M] <- Getcapcounts(y.true=y.true[1:M,1:J,1:K]) #intermediate object to derive n
  n <- Getncap(capcounts=capcounts[1:M],ID=ID[1:n.samples],y.state=y.state[1:M,1:J,1:K]) #number of captured individuals
  N <- sum(z[1:M])
})# end model
