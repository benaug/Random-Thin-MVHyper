#------------------------------------------------------------------
# Function for calculation detection rate
#------------------------------------------------------------------
GetKern <- nimbleFunction(
  run = function(s = double(1), sigma=double(0), 
                 X=double(2), J=double(0), z=double(0)){ 
    returnType(double(1))
    if(z==0) return(rep(0,J))
    if(z==1){
     d2 <- ((s[1]-X[1:J,1])^2 + (s[2]-X[1:J,2])^2)
     ans <- exp(-d2/(2*sigma^2))
     return(ans)
    }
  }
)
GetLam <- nimbleFunction(
  run = function(kern=double(1),lam0=double(0),J=double(0), z=double(0)){ 
    returnType(double(1))
    if(z==0) return(rep(0,J))
    if(z==1){
      ans <- lam0*kern
      return(ans)
    }
  }
)

dPoisson2D <- nimbleFunction(
  run = function(x = double(2), lam.p = double(1), lam.c=double(1), y.state = double(2), K2D = double(2), z = double(0),
                 log = integer(0)) {
    returnType(double(0))
    if(z==0){
      if(sum(x)>0){ #need this so z is not turned off if samples allocated to individual
        return(-Inf)
      }else{
        return(0)
      }
    }else{
      J <- nimDim(K2D)[1]
      K <- nimDim(K2D)[2]
      logProb <- matrix(NA,J,K)
      for(j in 1:J){
        for(k in 1:K){
          if(K2D[j,k]>0){
            if(y.state[j,k]==0){
              logProb[j,k] <- dpois(x[j,k], lambda = lam.p[j], log = TRUE)
            }else{
              logProb[j,k] <- dpois(x[j,k], lambda = lam.c[j], log = TRUE)
            }
          }
        }
      }
      return(sum(logProb))
    }
  }
)
#make dummy random vector generator to make nimble happy
rPoisson2D <- nimbleFunction(
  run = function(n = integer(0),lam.p = double(1), lam.c=double(1), y.state = double(2), K2D=double(2), z = double(0)) {
    returnType(double(2))
    J <- nimDim(lam.p)[1]
    K <- nimDim(lam.p)[2]
    out <- matrix(nrow=J,ncol=K,value=0)
    return(out)
  }
)

Getcapcounts <- nimbleFunction(
  run = function(y.true=double(3)){
    returnType(double(1))
    M <- nimDim(y.true)[1]
    J <- nimDim(y.true)[2]
    K <- nimDim(y.true)[3]
    capcounts=numeric(M, value = 0)
    for(i in 1:M){
      capcounts[i]=sum(y.true[i,1:J,1:K])
    }
    return(capcounts)
  }
)
Getncap <- nimbleFunction(
  run = function(capcounts=double(1),ID=double(1),y.state=double(3)){ #don't need ID, but nimble requires is it used in a function 
    returnType(double(0))
    M <- nimDim(capcounts)[1]
    nstate <- numeric(M, value = 0)
    for(i in 1:M){
      if(capcounts[i]>0){
        nstate[i] <- 1
      }
    }
    n.cap <- sum(nstate)
    return(n.cap)
  }
)

NimChoose <- nimbleFunction(
  run = function(n=double(0),k=double(0)){
    returnType(double(0))
    out <- (lfactorial(n)-lfactorial(k)-lfactorial(n-k))
    return(exp(out))
  }
)

#------------------------------------------------------------------
# Customer sampler to update latent IDs, and associated arrays
#------------------------------------------------------------------
IDSampler <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    M<-control$M
    J <- control$J
    K <- control$K
    K2D <- control$K2D
    n.samples <- control$n.samples
    this.j <- control$this.j
    this.k <- control$this.k
    y.ID <- control$y.ID
    n.jk <- control$n.jk
    denom.choose <- control$denom.choose
    calcNodes <- model$getDependencies(target)
  },
  run = function() {
    z <- model$z
    y.true <- model$y.true
    y.state <- model$y.state
    
    #precalculate match. Does sample l match individual i?
    match <- matrix(TRUE,nrow=n.samples,ncol=M) #start with all TRUE
    for(i in 1:M){#can match any individual
      if(z[i]==0){#unless z is off
        match[1:n.samples,i] <- FALSE
      }
    }
    
    #precalculate log likelihoods.
    ll.y <- array(0,dim=c(M,J,K))
    ll.y.thin <- matrix(0,J,K)
    for(i in 1:M){
      if(z[i]==1){
        for(j in 1:J){
          for(k in 1:K){
            if(K2D[j,k]>0){
              if(y.state[i,j,k]==0){
                ll.y[i,j,k] <- dpois(y.true[i,j,k],model$lam.p[i,j],log=TRUE)
              }else{
                ll.y[i,j,k] <- dpois(y.true[i,j,k],model$lam.c[i,j],log=TRUE)
              }
            }
          }
        }
      }
    }
    
    for(j in 1:J){
      for(k in 1:K){
        if(n.jk[j,k]>0){
          prod.choose <- 1
          for(i in 1:M){
            # if(z[i]==1){
            if(model$capcounts[i]>0){ #can screen out more 0's using capcounts here. Can't use for proposed likelihood below bc not updated
              prod.choose <- prod.choose*NimChoose(y.true[i,j,k],y.ID[i,j,k])
            }
          }
          ll.y.thin[j,k] <- log(prod.choose/denom.choose[j,k])
        }
      }
    }

    ll.y.cand <- ll.y
    ll.y.thin.cand <- ll.y.thin
    ID.curr <- model$ID

    ###update IDs
    ID.cand <- ID.curr
    y.true.cand <- y.true
    y.state.cand <- y.state
    for(l in 1:n.samples){#for all samples without known IDs
      propprobs <- model$lam.p[1:M,this.j[l]]
      for(i in 1:M){ #zero out nonmatches and z=0
        if(!match[l,i]){
          propprobs[i] <- 0
        }
      }
      denom <- sum(propprobs) #abort if propprobs sum to 0. No matches anywhere nearby.
      if(denom>0){
        propprobs <- propprobs/denom
        ID.cand[l] <- rcat(1,prob=propprobs)
        if(ID.cand[l]!=ID.curr[l]){
          swapped <- c(ID.curr[l],ID.cand[l])
          #new sample proposal probabilities
          forprob <- propprobs[swapped[2]]
          backprob <- propprobs[swapped[1]]
          #new y.true's - move sample from ID to ID.cand
          y.true.cand[ID.curr[l],this.j[l],this.k[l]] <- y.true[ID.curr[l],this.j[l],this.k[l]]-1
          y.true.cand[ID.cand[l],this.j[l],this.k[l]] <- y.true[ID.cand[l],this.j[l],this.k[l]]+1
          #update y.state
          #first guy
          for(k in 1:K){
            y.state.cand[swapped[1],this.j[l],k] <- 0
          }
          hasSamps <- sum(y.true.cand[swapped[1],this.j[l],])>0
          if(hasSamps){
            firstcap <- which(y.true.cand[swapped[1],this.j[l],]>0)[1]
            if(firstcap<K){
              for(k in (firstcap+1):K){
                y.state.cand[swapped[1],this.j[l],k] <- 1
              }
            }
          }
          #second guy
          for(k in 1:K){
            y.state.cand[swapped[2],this.j[l],k] <- 0
          }
          hasSamps <- sum(y.true.cand[swapped[2],this.j[l],])>0
          if(hasSamps){
            firstcap <- which(y.true.cand[swapped[2],this.j[l],]>0)[1]
            if(firstcap<K){
              for(k in (firstcap+1):K){
                y.state.cand[swapped[2],this.j[l],k] <- 1
              }
            }
          }
          
          #proposal likelihoods
          for(k in 1:K){#all k may change
            if(K2D[this.j[l],k]>0){
              if(y.state.cand[swapped[1],this.j[l],k]==0){
                ll.y.cand[swapped[1],this.j[l],k] <- dpois(y.true.cand[swapped[1],this.j[l],k],model$lam.p[swapped[1],this.j[l]],log=TRUE)
              }else{
                ll.y.cand[swapped[1],this.j[l],k] <- dpois(y.true.cand[swapped[1],this.j[l],k],model$lam.c[swapped[1],this.j[l]],log=TRUE)
              }
              if(y.state.cand[swapped[2],this.j[l],k]==0){
                ll.y.cand[swapped[2],this.j[l],k] <- dpois(y.true.cand[swapped[2],this.j[l],k],model$lam.p[swapped[2],this.j[l]],log=TRUE)
              }else{
                ll.y.cand[swapped[2],this.j[l],k] <- dpois(y.true.cand[swapped[2],this.j[l],k],model$lam.c[swapped[2],this.j[l]],log=TRUE)
              }
            }
          }
          prod.choose <- 1
          for(i in 1:M){
            if(z[i]==1){#screening out some 0's by skipping z=0 here.
              prod.choose <- prod.choose*NimChoose(y.true.cand[i,this.j[l],this.k[l]],y.ID[i,this.j[l],this.k[l]])
            }
          }
          ll.y.thin.cand[this.j[l],this.k[l]] <- log(prod.choose/denom.choose[this.j[l],this.k[l]])

          #select sample to move proposal probabilities
          #P(select a sample of this type (not ID'd) for this ID)*P(select this j-k|sample of this type and this ID)
          #n.samples cancels out in MH ratio. Including for clarity
          focalprob <- (sum(ID.curr==swapped[1])/n.samples)*
            (y.true[swapped[1],this.j[l],this.k[l]] - y.ID[swapped[1],this.j[l],this.k[l]])/sum(y.true[swapped[1],1:J,1:K] - y.ID[swapped[1],1:J,1:K])
          focalbackprob <- (sum(ID.cand==swapped[2])/n.samples)*
            (y.true.cand[swapped[2],this.j[l],this.k[l]]- y.ID[swapped[2],this.j[l],this.k[l]])/sum(y.true.cand[swapped[2],1:J,1:K] - y.ID[swapped[2],1:J,1:K])

          #sum log likelihoods and do MH step
          lp_initial <- sum(ll.y[swapped[1],this.j[l],])+
            sum(ll.y[swapped[2],this.j[l],])+
            ll.y.thin[this.j[l],this.k[l]]
          lp_proposed <- sum(ll.y.cand[swapped[1],this.j[l],])+
            sum(ll.y.cand[swapped[2],this.j[l],])+
            ll.y.thin.cand[this.j[l],this.k[l]]
          log_MH_ratio <- (lp_proposed+log(backprob)+log(focalbackprob)) - (lp_initial+log(forprob)+log(focalprob))
          accept <- decide(log_MH_ratio)

          if(accept){
            y.true[swapped[1],this.j[l],this.k[l]] <- y.true.cand[swapped[1],this.j[l],this.k[l]]
            y.true[swapped[2],this.j[l],this.k[l]] <- y.true.cand[swapped[2],this.j[l],this.k[l]]
            y.state[swapped[1],this.j[l],] <- y.state.cand[swapped[1],this.j[l],]
            y.state[swapped[2],this.j[l],] <- y.state.cand[swapped[2],this.j[l],]
            ll.y[swapped[1],this.j[l],] <- ll.y.cand[swapped[1],this.j[l],]
            ll.y[swapped[2],this.j[l],] <- ll.y.cand[swapped[2],this.j[l],]
            ll.y.thin[this.j[l],this.k[l]] <- ll.y.thin.cand[this.j[l],this.k[l]]
            ll.y.thin[this.j[l],this.k[l]] <- ll.y.thin.cand[this.j[l],this.k[l]]
            ID.curr[l] <- ID.cand[l]
          }else{
            #set these back.
            y.true.cand[swapped[1],this.j[l],this.k[l]] <- y.true[swapped[1],this.j[l],this.k[l]]
            y.true.cand[swapped[2],this.j[l],this.k[l]] <- y.true[swapped[2],this.j[l],this.k[l]]
            ID.cand[l] <- ID.curr[l]
          }
        }
      }
    }
    
    #put everything back into the model$stuff
    model$y.true <<- y.true
    model$y.state <<- y.state
    model$ID <<- ID.curr
    model.lp.proposed <- model$calculate(calcNodes) #update logprob
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list( reset = function () {} )
)