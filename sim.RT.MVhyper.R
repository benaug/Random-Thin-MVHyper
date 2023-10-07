e2dist <- function (x, y){
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

sim.RT.MVhyper <-
  function(N=NA,lam0=NA,sigma=NA,p0=NA,lambda=NA,theta.d=NA,K=10,X=NA,buff=NA,obstype="poisson",
           n.select=NA,K2D=NA){
    library(abind)
    
    # simulate a population of activity centers
    X <- as.matrix(X)
    xlim <- c(min(X[,1]),max(X[,1]))+c(-buff,buff)
    ylim <- c(min(X[,2]),max(X[,2]))+c(-buff,buff)
    s <- cbind(runif(N, xlim[1],xlim[2]), runif(N,ylim[1],ylim[2]))
    D <- e2dist(s,X)
    J <- nrow(X)
    
    #trap operation
    if(any(is.na(K2D))){
      K2D <- matrix(1,J,K)
    }
    
    # Capture individuals
    y.true <-array(0,dim=c(N,J,K))
    if(obstype=="bernoulli"){
      if(is.na(lam0))stop("must provide lam0 for bernoulli obstype")
      pd <- p0 * exp(-D * D/(2 * sigma * sigma))
      for(i in 1:N){
        for(j in 1:J){
          for(k in 1:K){
            y.true[i,j,k] <- rbinom(1,1,pd[i,j]*K2D[j,k])
          }
        }
      }
    }else if(obstype=="poisson"){
      if(is.na(lam0))stop("must provide lam0 for poisson obstype")
      lamd <- lam0*exp(-D*D/(2*sigma*sigma))
      for(i in 1:N){
        for(j in 1:J){
          for(k in 1:K){
            y.true[i,j,k] <- rpois(1,lamd[i,j]*K2D[j,k])
          }
        }
      }
    }else if(obstype=="negbin"){
      if(is.na(lam0))stop("must provide lam0 for negbin obstype")
      if(is.na(theta.d))stop("Must provide theta.d for negbin obstype")
      lamd <- lam0*exp(-D*D/(2*sigma*sigma))
      for(i in 1:N){
        for(j in 1:J){
          for(k in 1:K){
            if(K2D[j,k]>0){
              y.true[i,j,k] <- rnbinom(1,mu=lamd[i,j],size=theta.d)
            }
          }
        }
      } 
    }else if(obstype=="hurdleZTPois"){
      if(is.na(p0))stop("must provide p0 for siteUse_ZTpois obstype")
      if(is.na(lambda))stop("must provide lambda for siteUse_ZTpois obstype")
      library(VGAM)
      pd <- p0*exp(-D*D/(2*sigma*sigma))
      y.det <- array(0,dim=c(N,J,K))
      for(i in 1:N){
        for(j in 1:J){
          for(k in 1:K){
            y.det[i,j,k] <- rbinom(1,1,pd[i,j])
            if(y.det[i,j,k]==1){
              y.true[i,j,k] <- rzapois(1,lambda,pobs0=0)
            }
          }
        }
      }
    }else{
      stop("obstype not recognized")
    }
    
    y.true.full <- y.true
    n.cap <- sum(rowSums(y.true.full)>0)
    
    captured <- which(rowSums(y.true)>0)
    y.true <- y.true[captured,,]
    n.cap <- nrow(y.true)
   
    #split sightings into ID'd or not ID'd
    n.samples <- sum(y.true)
    ID <- rep(NA,n.samples)
    y.ID <- y.true
    y.noID <- array(0,dim=c(n.samples,J,K))
    idx <- 1
    
    #loop through traps and occasions subsampling each
    for(j in 1:J){ #then traps
      for(k in 1:K){ #then occasions
        if(sum(y.true[,j,k])>0){ #is there at least one sample here?
          these.guys <- which(y.true[,j,k]>0)
          these.samps <- y.true[these.guys,j,k]
          ids <- rep(these.guys,times=these.samps)
          if(length(ids)<=n.select){#if fewer than or the same number of samples as n.select, keep them all
            pick <- 1:length(ids)
          }else{#otherwise, subsample
            pick <- sample(1:length(ids),size=n.select,replace=FALSE)
          }
          remaining.ids=ids[-pick]
          if(length(remaining.ids)>0){
            for(i in remaining.ids){
              y.noID[idx,j,k] <- 1 #add to not identified data
              y.ID[i,j,k] <- y.ID[i,j,k]-1 #subtract from identified data
              ID[idx] <- i
              idx <- idx+1
            }
          }
        }
      }
    }
    
    n.samples.noID <- sum(rowSums(y.noID))
    y.noID <- y.noID[1:n.samples.noID,,]
    ID <- ID[1:n.samples.noID]
    
    #reassemble to check for correctness
    y.test <- y.ID
    for(i in 1:n.samples.noID){
      y.test[ID[i],,] <- y.test[ID[i],,]+y.noID[i,,]
    }
    
   if(!all(y.test==y.true))stop("Error reassembling data. Bug in simulator.") #shouldn't happen 

    #reorder ID's if any y.ID inds end up with 0 samples after subsampling.
    #only needed to compare posterior ID match probs to correct ID #
    allIDs <- 1:nrow(y.ID) #all guys in y.ID with samples before subsampling
    ID.map <- allIDs #what we'll use to map old IDs to new IDs
    fix <- which(rowSums(y.ID)==0) #guys subsampled out
    if(length(fix)>0){
      ID.map[fix] <- (max(allIDs)+1):(max(allIDs)+length(fix)) #give these guys new ID numbers not already used
      #then subtract the removed guys from ID numbers of other guy's samples
      for(i in 1:length(fix)){
        ID.map[ID.map>=fix[i]] <- ID.map[ID.map>=fix[i]]-1
        fix <- fix-1
      }
      ID <- ID.map[ID]
    }
    
    #remove inds subsampled out of y.ID
    IDd <- which(rowSums(y.ID)!=0)
    n.ID <- length(IDd)
    y.ID <- y.ID[IDd,,]
    
    #observed capture data can be represented by site and occasion of each count member
    this.j <- this.k <- rep(NA,n.samples.noID)
    for(i in 1:n.samples.noID){
      tmp <- which(y.noID[i,,]==1,arr.ind=TRUE)
      this.j[i] <- tmp[1]
      this.k[i] <- tmp[2]
    }
    
    #reorder ID, this.j, this.k by ID (not required)
    ord <- order(ID)
    ID <- ID[ord]
    this.j <- this.j[ord]
    this.k <- this.k[ord]
    

    out <- list(y.true=y.true,y.ID=y.ID,this.j=this.j,this.k=this.k,
              X=X,K=K,K2D=K2D,buff=buff,s=s,xlim=xlim,ylim=ylim,
              ID=ID,n.ID=n.ID,y.true.full=y.true.full,n.cap=n.cap)
    return(out)
  }