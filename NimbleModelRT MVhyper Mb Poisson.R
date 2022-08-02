NimModel <- nimbleCode({
  #detection function priors
  lam0.p~dunif(0,20)
  lam0.c~dunif(0,20)
  sigma~dunif(0,20)
  #data augmentation prior
  psi~dunif(0,1)
  #likelihoods (except for s priors)
  for(i in 1:M){
    z[i] ~ dbern(psi)
    s[i,1] ~ dunif(xlim[1],xlim[2])
    s[i,2] ~ dunif(ylim[1],ylim[2])
    kern[i,1:J] <- GetKern(s = s[i,1:2], X = X[1:J,1:2], J=J,sigma=sigma, z=z[i])
    lam.p[i,1:J] <- GetLam(kern = kern[i,1:J],lam0=lam0.p,J=J,z=z[i])
    lam.c[i,1:J] <- GetLam(kern = kern[i,1:J],lam0=lam0.c,J=J,z=z[i])
    y.true[i,1:J,1:K] ~ dPoisson2D(lam.p[i,1:J],lam.c[i,1:J],y.state[i,1:J,1:K],
                                   K2D[1:J,1:K],z=z[i])  # Model for complete capture histories
  }
  #derived variables
  capcounts[1:M] <- Getcapcounts(y.true=y.true[1:M,1:J,1:K]) #intermediate object to derive n
  n <- Getncap(capcounts=capcounts[1:M],ID=ID[1:n.samples],y.state=y.state[1:M,1:J,1:K]) #number of captured individuals
  N <- sum(z[1:M])
})# end model
