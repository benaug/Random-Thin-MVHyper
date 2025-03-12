e2dist<- function (x, y){
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

init.RT.MVhyper.Mb=function(data,inits=NA,M=NA,obstype="poisson"){
  library(abind)
  y.ID <- data$y.ID
  this.j <- data$this.j
  this.k <- data$this.k
  X <- as.matrix(data$X)
  J <- nrow(X)
  K <- dim(y.ID)[3]
  buff <- data$buff
  n.ID <- data$n.ID
  K2D <- data$K2D
  
  #data checks
  if(length(dim(y.ID))!=3){
    stop("dim(y.ID) must be 3. Reduced to 2 during initialization")
  }

  buff <- data$buff
  xlim <- c(min(X[,1]),max(X[,1]))+c(-buff, buff)
  ylim <- c(min(X[,2]),max(X[,2]))+c(-buff, buff)
  
  ##pull out initial values
  lam0.p <- inits$lam0.p
  lam0.c <- inits$lam0.c
  sigma <- inits$sigma
  
  #initialize IDs
  #This initialization algorithm matches samples with consistent G.noID to same ID if they are caught
  #at same trap. Requires larger data augmentation than necessary for sampling. Improvement to implement:
  #combine if consistent G.noID and caught at traps within a distance consistent with initial sigma.
  n.samples=length(this.j)
  ID=rep(NA,n.samples)
  nextID=n.ID+1
  
  library(abind)
  y.ID=abind(y.ID,array(0,dim=c(M-n.ID,J,K)),along=1)
  y.true=y.ID
  for(l in 1:n.samples){
    sametrap=y.true[,this.j[l],this.k[l]]>0
    if(any(sametrap)){
      ID[l]=which(sametrap)[1]
      y.true[ID[l],this.j[l],this.k[l]]=y.true[ID[l],this.j[l],this.k[l]]+1
    }else{#must be new ID
      ID[l]=nextID
      y.true[ID[l],this.j[l]]=y.true[ID[l],this.j[l],this.k[l]]+1
      nextID=nextID+1
    }
    if(nextID>M)stop("Need to raise M to initialize data.")
  }
  
  #initialize z
  z=1*(rowSums(y.true)>0)

  #intialize s
  y.true2D=apply(y.true,c(1,2),sum)
  s<- cbind(runif(M,xlim[1],xlim[2]), runif(M,ylim[1],ylim[2])) #assign random locations
  idx=which(rowSums(y.true2D)>0) #switch for those actually caught
  for(i in idx){
    trps<- matrix(X[y.true2D[i,]>0,1:2],ncol=2,byrow=FALSE)
    if(nrow(trps)>1){
      s[i,]<- c(mean(trps[,1]),mean(trps[,2]))
    }else{
      s[i,]<- trps
    }
  }
  
  #get state matrix
  y.state=array(0,dim=c(M,J,K))
  for(i in 1:M){
    for(j in 1:J){
      if(sum(y.true[i,j,])>0){
        firstcap=which(y.true[i,j,]>0)[1]
        if(firstcap<K){
          y.state[i,j,(firstcap+1):K]=1
        }
      }
    }
  }
  
  D=e2dist(s, X)
  lam.p=lam0.p * exp(-D * D/(2 * sigma * sigma))
  lam.c=lam0.c * exp(-D * D/(2 * sigma * sigma))
  ll.y=y.true*0
  if(obstype=="poisson"){
    for(i in 1:M){
      if(z[i]==1){
        for(j in 1:J){
          for(k in 1:K){
            if(y.state[i,j,k]==0){
              ll.y[i,j,k]=dpois(y.true[i,j,k],K2D[j,k]*lam.p[i,j],log=TRUE)
            }else{
              ll.y[i,j,k]=dpois(y.true[i,j,k],K2D[j,k]*lam.c[i,j],log=TRUE)
            }
          }
        }
      }
    }
  }else if(obstype=="negbin"){
    theta.d=inits$theta.d
    ll.y=y.true*0
    for(i in 1:M){
      if(z[i]==1){
        for(j in 1:J){
          for(k in 1:K){
            if(y.state[i,j,k]==0){
              ll.y[i,j,k]=dnbinom(y.true[i,,k],mu=lam.p[i,j],size=theta.d*K2D[j,k],log=TRUE)
            }else{
              ll.y[i,j,k]=dnbinom(y.true[i,,k],mu=lam.c[i,j],size=theta.d*K2D[j,k],log=TRUE)
            }
          }
        }
      }
    }
  }else if(obstype=="hurdleZTPois"){
    p0.p=inits$p0.p
    p0.c=inits$p0.c
    lambda=inits$lambda
    pd.p <-p0.p * exp(-D * D/(2 * sigma * sigma))
    pd.c <-p0.c * exp(-D * D/(2 * sigma * sigma))
    for(i in 1:M){
      if(z[i]==1){
        for(j in 1:J){
          for(k in 1:K){
            if(K2D[j,k]>0){
              if(y.state[i,j,k]==0){
                if(y.true[i,j,k]==0){
                  ll.y[i,j,k]=log(1-pd.p[i,j])
                }else{
                  ll.y[i,j,k]=log(pd.p[i,j]) + log(dpois(y.true[i,j,k],lambda=lambda)/(1-exp(-lambda)))
                }
              }else{
                if(y.true[i,j,k]==0){
                  ll.y[i,j,k]=log(1-pd.c[i,j])
                }else{
                  ll.y[i,j,k]=log(pd.c[i,j]) + log(dpois(y.true[i,j,k],lambda=lambda)/(1-exp(-lambda)))
                }
              }
            }
          }
        }
      }
    }
  }else{
    stop("obstype not recognized")
  }
  
  if(!is.finite(sum(ll.y)))stop("Starting observation model likelihood not finite. Possible error in K2D (if supplied by user) or problem initializing data.")
  
  return(list(s=s,z=z,ID=ID,y.ID=y.ID,y.true=y.true,y.state=y.state,K2D=K2D,
         n.samples=n.samples,this.j=this.j,this.k=this.k,xlim=xlim,ylim=ylim))

}