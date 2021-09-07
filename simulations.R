library(parallel)
library(SQUAREM)

source("/home/emil/projects/EMasso/scripts/EMasso/regression_from_design.R")
source("/home/emil/projects/EMasso/scripts/EMasso/EMassoR.R")
source("/home/emil/projects/EMasso/scripts/EMasso/EMassoSimR.R")

args<-commandArgs(trailingOnly = T)

Nsimulations<-as.numeric(args[1])

emFreqAdmix <- function( x, GL, anc ){
  f <- x
  tol <- 1e-4
  anc[anc<tol] <- tol
  anc[anc > 1-tol] <- 1-tol
  Q <- rbind( anc, 1-anc )
  
  aNorm <- t( Q ) %*% f
  bNorm <- t( Q ) %*% ( 1-f )
  
  Htop<-GL[,2]*2*(aNorm*(1-aNorm))+GL[,3]*2*(aNorm)**2
  Hbottom<-GL[,1]*(1-aNorm)**2+GL[,2]*2*(aNorm*(1-aNorm))+GL[,3]*(aNorm)**2
  H<-Htop/Hbottom
  
  updateF <- function( k ){
    ag <- H*Q[k, ]*f[k]/aNorm
    bg <- ( 2-H )*Q[k, ]*( 1-f[k] )/bNorm
    fnew <- sum( ag )/( sum( ag ) + sum( bg ) )
    fnew
  }
  f <- sapply( 1:2, updateF )
  
  f[f<tol] <- tol
  f[1-f > 1-tol] <- 1-tol
  f
}

emFreqAdmixGeno <- function( x, geno, a ){
  f <- x
  q <- rbind( a, 1-a )
  tol <- 1e-4
  aNorm <- t( q ) %*% f
  bNorm <- t( q ) %*% ( 1-f )
  updateF <- function( k ){
    ag <- geno*q[k, ]*f[k]/aNorm
    bg <- ( 2-geno )*q[k, ]*( 1-f[k] )/bNorm
    fnew <- sum( ag )/( sum( ag ) + sum( bg ) )
    fnew
  }
  f <- sapply( 1:2, updateF )
  
  f[f<tol] <- tol
  f[1-f > 1-tol] <- 1-tol
  f
}


getLikes<-function(x,d=2,e=0.01,norm=FALSE){
  n<-length(x)
  dep<-rpois(n,d)
  nA<-rbinom(n,dep,c(e,0.5,1-e)[x+1])
  res<-rbind(dbinom(nA,dep,e),dbinom(nA,dep,0.5),dbinom(nA,dep,1-e))
  if(norm){
    res<-t(t(res)/colSums(res))
  }
  res[is.na(res)]<-1/3
  
  return(list(res=res,depth=dep))
}

## estimate freq from genotype likelihoods
estPem<-function(like){ #numeric optimazition EM
  pk<-0.001 #start
  for(tal in 1:20){
    w0<-like[,1]*(1-pk)^2
    w1<-like[,2]*2*pk*(1-pk)
    w2<-like[,3]*pk^2
    pk<-mean((w1+2*w2)/(2*(w0+w1+w2)))
  }
  pk
}


qqpemil2<-function(x,ci=TRUE,add=FALSE,ylab="Observed log10(p-value)",xlab="Expected log10(p-value)",maxLogP,...){
  x<-x[!is.na(x)]
  if(!missing(maxLogP)){
    x[x<10^-maxLogP]<-10^-maxLogP
  }
  N<-length(x)
  chi1<-qchisq(1-x,1)
  x<-sort(x)
  lambda<-round(median(chi1)/qchisq(0.5,1),2)
  e<- -log((1:N-0.5)/N,10)
  if(add){
    points(e,-log(x,10),pch=20,...)
  } else{
    plot(e,-log(x,10),ylab=ylab,xlab=xlab,pch=20,...)
    abline(0,1,col=2,lwd=2)
  }
    
  if(ci){
    c95<-qbeta(0.95,1:N,N-(1:N)+1)
    c05<-qbeta(0.05,1:N,N-(1:N)+1)
    lines(e,-log(c95,10))
    lines(e,-log(c05,10))
  }
}


fun <- function(seed,delta,gamma=0,beta=0,keepInformative=FALSE,estimateMaf=0,par){
  
  g <- rbinom(par$Nind,2,par$fInd)  
  y <- rnorm(par$Nind,sd=par$sd)+beta*g + gamma*par$Q[1,]  
  depthGroup<-c(rep(par$depthGroup[1],par$Nind/2),rep(par$depthGroup[2],par$Nind/2))  
  prob<-1/(1+exp(-delta*y))  
  dep <- ifelse(rbinom(par$Nind,1,prob)==1,par$depthGroup[1],par$depthGroup[2])  
  likes<-getLikes(g,d=dep,norm=T,e=0.01)
  GL <- likes$res  
  keep<-likes$depth>0
  
  if(!(estimateMaf%in%c(0,1,2))){
    stop("estimateMaf has to be 0, 1 or 2")    
  }
  
  if(keepInformative){
    
    if(estimateMaf==0){
      maf <- squarem( c( 0.1, 0.1 ), emFreqAdmix, GL = t(GL[,keep]), anc = par$Q[1,keep] )$par
      f<-estPem(t(GL[,keep]))
    } else if(estimateMaf==1){
      maf<- squarem( c( 0.1, 0.1 ), emFreqAdmixGeno, geno = g[keep], a = par$Q[1,keep] )$par
      f <- sum(g[keep])/sum(keep)/2
    } else if(estimateMaf==2){
      maf<-par$freq
      f<-mean(par$Q[1,keep]*maf[1] + ( 1-par$Q[1,keep] )*maf[2])
    }
    
    fInd <- par$Q[1,keep]*maf[1] + ( 1-par$Q[1,keep] )*maf[2]
    
    PPf <- GL[,keep]*c((1-f)^2,2*f*(1-f),f^2)
    PPf <- PPf/rep(colSums(PPf),each=3)
    Dosf<- colSums(PPf*c(0,1,2))
    
    PPind <- GL[,keep]*rbind((1-fInd)^2,2*fInd*(1-fInd),fInd^2)
    PPind <- PPind/rep(colSums(PPind),each=3)
    Dosind<- colSums(PPind*c(0,1,2))
    
    g<-g[keep]
    y<-y[keep]
    dep<-dep[keep]
    Q<-par$Q[,keep]
    depthGroup<-depthGroup[keep]
    
  } else{
    
    Q<-par$Q
    
    if(estimateMaf==0){
      maf <- squarem( c( 0.1, 0.1 ), emFreqAdmix, GL = t(GL), anc = Q[1,] )$par
      f<-estPem(t(GL))
    } else if(estimateMaf==1){
      maf<- squarem( c( 0.1, 0.1 ), emFreqAdmixGeno, geno = g, a = Q[1,] )$par
      f <- sum(g)/par$Nind/2
    } else if(estimateMaf==2){
      maf<-par$freq
      f<-mean(Q[1,]*maf[1] + ( 1-Q[1,] )*maf[2])
    }
    
    fInd <- Q[1,]*maf[1] + ( 1-Q[1,] )*maf[2]
    
    PPf <- GL*c((1-f)^2,2*f*(1-f),f^2)
    PPf <- PPf/rep(colSums(PPf),each=3)
    Dosf<- colSums(PPf*c(0,1,2))
    
    PPind <- GL*rbind((1-fInd)^2,2*fInd*(1-fInd),fInd^2)
    PPind <- PPind/rep(colSums(PPind),each=3)
    Dosind<- colSums(PPind*c(0,1,2))
    
  }  
  
  fadmix<-tapply(g,Q[1,],function(x) mean(x)/2)
  fdepth<-tapply(g,depthGroup,function(x) mean(x)/2)
    
  if(all(Q[1,]==1) | all(Q[1,]==0)){
    list(
      f=f,
      fadmix=fadmix,
      fdepth=fdepth,
      PPf=EMasso(GP=t(PPf),ys=y),
      PPind=EMasso(GP=t(PPind),ys=y),
      Dosf=summary(lm(y~Dosf)),
      Dosind=summary(lm(y~Dosind)),
      geno=summary(lm(y~g))
    )
    
  } else{
    list(
      f=f,
      fadmix=fadmix,
      fdepth=fdepth,
      PPf=EMasso(GP=t(PPf),ys=y,covs=cbind(1,Q[1,])),
      PPind=EMasso(GP=t(PPind),ys=y,covs=cbind(1,Q[1,])),
      Dosf=summary(lm(y~Dosf+Q[1,])),
      Dosind=summary(lm(y~Dosind+Q[1,])),
      geno=summary(lm(y~g+Q[1,]))
    )
  }
  
}


fun2 <- function(delta,gamma=0,beta=0,keepInformative=FALSE,estimateMaf=0,par){
  cat("delta",delta," gamma",gamma," beta",beta,"\n")
  res <- parallel::mclapply(1:par$simulations, function(x) fun(x,delta=delta,gamma=gamma,beta=beta,keepInformative = keepInformative,estimateMaf=estimateMaf,par=par),mc.cores=10)
  
  if(all(par$Q[1,]==1) | all(par$Q[1,]==0)){
    
    powInd <- mean( sapply(res,function(x) x$PPind$test[3] ) < par$threshold )
    powDosInd <- mean( sapply(res,function(x) x$Dosind$coefficients[8] ) < par$threshold )
    
    powF <- mean( sapply(res,function(x) x$PPf$test[3] ) < par$threshold )
    powDosF <- mean( sapply(res,function(x) x$Dosf$coefficients[8] ) < par$threshold )
    
    powG <- mean( sapply(res,function(x) x$geno$coefficients[8] ) < par$threshold )
    
    diffF<-mean( sapply(res,function(x) x$PPf$full$par[1]-beta ) )
    diffInd<- mean( sapply(res,function(x) x$PPind$full$par[1]-beta ) )
    
    diffDosF<-mean( sapply(res,function(x) x$Dosf$coefficients[2]-beta ) )
    diffDosInd<- mean( sapply(res,function(x) x$Dosind$coefficients[2]-beta ) )
    
    diffGeno<- mean( sapply(res,function(x) x$geno$coefficients[2]-beta ) )
    
    betaF<-mean(sapply(res,function(x) x$PPf$full$par[1]))
    betaDosF<-mean(sapply(res,function(x) x$Dosf$coefficients[2]))
    betaGeno<-mean(sapply(res,function(x) x$geno$coefficients[2]))
    
    alphaF<-mean(sapply(res,function(x) x$PPf$full$par[2]))
    alphaDosF<-mean(sapply(res,function(x) x$Dosf$coefficients[1]))
    alphaGeno<-mean(sapply(res,function(x) x$geno$coefficients[1]))
    
    f<- sapply(res,function(x) x$f )
    
    fadmix<- sapply(res,function(x) x$fadmix )
    fdepth<- sapply(res,function(x) x$fdepth )
    
  } else{
    
    powInd <- mean( sapply(res,function(x) x$PPind$test[3] ) < par$threshold )
    powDosInd <- mean( sapply(res,function(x) x$Dosind$coefficients[11] ) < par$threshold )
    
    powF <- mean( sapply(res,function(x) x$PPf$test[3] ) < par$threshold )
    powDosF <- mean( sapply(res,function(x) x$Dosf$coefficients[11] ) < par$threshold )
    
    powG <- mean( sapply(res,function(x) x$geno$coefficients[11] ) < par$threshold )
    
    diffF<-mean( sapply(res,function(x) x$PPf$full$par[1]-beta ) )
    diffInd<- mean( sapply(res,function(x) x$PPind$full$par[1]-beta ) )
    
    diffDosF<-mean( sapply(res,function(x) x$Dosf$coefficients[2]-beta ) )
    diffDosInd<- mean( sapply(res,function(x) x$Dosind$coefficients[2]-beta ) )
    
    diffGeno<- mean( sapply(res,function(x) x$geno$coefficients[2]-beta ) )
    
    betaF<-mean(sapply(res,function(x) x$PPf$full$par[1]))
    betaDosF<-mean(sapply(res,function(x) x$Dosf$coefficients[2]))
    betaGeno<-mean(sapply(res,function(x) x$geno$coefficients[2]))
    
    alphaF<-mean(sapply(res,function(x) x$PPf$full$par[2]))
    alphaDosF<-mean(sapply(res,function(x) x$Dosf$coefficients[1]))
    alphaGeno<-mean(sapply(res,function(x) x$geno$coefficients[1]))
    
    f<- sapply(res,function(x) x$f )
    
    fadmix<- sapply(res,function(x) x$fadmix )
    fdepth<- sapply(res,function(x) x$fdepth )
    
  }
  
  return(list(powInd=powInd,powDosInd=powDosInd,powF=powF,powDosF=powDosF,powG=powG,
              diffF=diffF,diffInd=diffInd,diffDosF=diffDosF,diffDosInd=diffDosInd,
              diffGeno=diffGeno,betaF=betaF,betaDosF=betaDosF,betaGeno=betaGeno,alphaF=alphaF,alphaDosF=alphaDosF,alphaGeno=alphaGeno,f=f,
              fadmix=fadmix,fdepth=fdepth))
  
}


funBin <- function(seed,delta,gamma=0,beta=0,keepInformative=FALSE,estimateMaf=0,par){
  
  g <- rbinom(par$Nind,2,par$fInd)  
  yBin <- rbinom(n = par$Nind, size = 1, prob = plogis(rnorm(par$Nind)+beta*g + gamma*par$Q[1,] ))  
  dep <- numeric(par$Nind)
  dep[yBin==0]<-1
  dep[yBin==1]<-4  
  likes<-getLikes(g,d=dep,norm=T,e=0.01)
  GL <- likes$res  
  keep<-likes$depth>0
  
  if(!(estimateMaf%in%c(0,1,2))){
    stop("estimateMaf has to be 0, 1 or 2")    
  }
  
  if(keepInformative){
    
    if(estimateMaf==0){
      maf <- squarem( c( 0.1, 0.1 ), emFreqAdmix, GL = t(GL[,keep]), anc = par$Q[1,keep] )$par
      f<-estPem(t(GL[,keep]))
    } else if(estimateMaf==1){
      maf<-par$freq
      f <- sum(g[keep])/sum(keep)/2
    } else if(estimateMaf==2){
      maf<-par$freq
      f<-mean(par$Q[1,keep]*maf[1] + ( 1-par$Q[1,keep] )*maf[2])
    }
    
    fInd <- par$Q[1,keep]*maf[1] + ( 1-par$Q[1,keep] )*maf[2]
    
    PPf <- GL[,keep]*c((1-f)^2,2*f*(1-f),f^2)
    PPf <- PPf/rep(colSums(PPf),each=3)
    Dosf<- colSums(PPf*c(0,1,2))
    
    PPind <- GL[,keep]*rbind((1-fInd)^2,2*fInd*(1-fInd),fInd^2)
    PPind <- PPind/rep(colSums(PPind),each=3)
    Dosind<- colSums(PPind*c(0,1,2))
    
    g<-g[keep]
    yBin<-yBin[keep]
    dep<-dep[keep]
    Q<-par$Q[,keep]
    
    
  } else{
    
    Q<-par$Q
    
    if(estimateMaf==0){
      maf <- squarem( c( 0.1, 0.1 ), emFreqAdmix, GL = t(GL), anc = Q[1,] )$par
      f<-estPem(t(GL))
    } else if(estimateMaf==1){
      maf<-par$freq
      f <- sum(g)/par$Nind/2
    } else if(estimateMaf==2){
      maf<-par$freq
      f<-mean(Q[1,]*maf[1] + ( 1-Q[1,] )*maf[2])
    }
    
    fInd <- Q[1,]*maf[1] + ( 1-Q[1,] )*maf[2]
    
    PPf <- GL*c((1-f)^2,2*f*(1-f),f^2)
    PPf <- PPf/rep(colSums(PPf),each=3)
    Dosf<- colSums(PPf*c(0,1,2))
    
    PPind <- GL*rbind((1-fInd)^2,2*fInd*(1-fInd),fInd^2)
    PPind <- PPind/rep(colSums(PPind),each=3)
    Dosind<- colSums(PPind*c(0,1,2))
    
  }
  
  if(all(Q[1,]==1) | all(Q[1,]==0)){
    list(
      f=f,
      PPf=EMasso(GP=t(PPf),ys=yBin,qu=F),
      PPind = EMasso(GP=t(PPind),ys=yBin,qu=F),
      Dosf=summary(glm(yBin~Dosf,family="binomial")),
      Dosind=summary(glm(yBin~Dosind,family="binomial")),
      geno=summary(glm(yBin~g,family="binomial"))
    )
  } else{
    list(
      f=f,
      PPf=EMasso(GP=t(PPf),ys=yBin,covs=cbind(1,Q[1,]),qu=F),
      PPind = EMasso(GP=t(PPind),ys=yBin,covs=cbind(1,Q[1,]),qu=F),
      Dosf=summary(glm(yBin~Dosf+Q[1,],family="binomial")),
      Dosind=summary(glm(yBin~Dosind+Q[1,],family="binomial")),
      geno=summary(glm(yBin~g+Q[1,],family="binomial"))
    )
  }
}


funBin2 <- function(delta,gamma=0,beta=0,keepInformative=FALSE,estimateMaf=0,par){
  cat("delta",delta," gamma",gamma," beta",beta,"\n")
  res <- parallel::mclapply(1:par$simulations, function(x) funBin(x,delta=delta,gamma=gamma,beta=beta,keepInformative = keepInformative,estimateMaf=estimateMaf,par=par),mc.cores=10)
  
  if(all(par$Q[1,]==1) | all(par$Q[1,]==0)){
    
    powInd <- mean( sapply(res,function(x) x$PPind$test[3] ) < par$threshold )
    powDosInd <- mean( sapply(res,function(x) x$Dosind$coefficients[8] ) < par$threshold )
    
    
    powF <- mean( sapply(res,function(x) x$PPf$test[3] ) < par$threshold )
    powDosF <- mean( sapply(res,function(x) x$Dosf$coefficients[8] ) < par$threshold )
    
    
    powG <- mean( sapply(res,function(x) x$geno$coefficients[8] ) < par$threshold )
    
    diffF<-mean( sapply(res,function(x) exp(x$PPf$full$par[1])-exp(beta) ) )
    diffInd<- mean( sapply(res,function(x) exp(x$PPind$full$par[1])-exp(beta) ) )
    
    diffDosF<-mean( sapply(res,function(x) exp(x$Dosf$coefficients[2])-exp(beta) ) )
    diffDosInd<- mean( sapply(res,function(x) exp(x$Dosind$coefficients[2])-exp(beta) ) )
    
    diffGeno<-mean( sapply(res,function(x) exp(x$geno$coefficients[2])-exp(beta) ) )
    
    betaF<-mean(sapply(res,function(x) x$PPf$full$par[1]))
    betaDosF<-mean(sapply(res,function(x) x$Dosf$coefficients[2]))
    betaGeno<-mean(sapply(res,function(x) x$geno$coefficients[2]))
    
    alphaF<-mean(sapply(res,function(x) x$PPf$full$par[2]))
    alphaDosF<-mean(sapply(res,function(x) x$Dosf$coefficients[1]))
    alphaGeno<-mean(sapply(res,function(x) x$geno$coefficients[1]))
    
    f<- sapply(res,function(x) x$f )
    
  } else{
    
    powInd <- mean( sapply(res,function(x) x$PPind$test[3] ) < par$threshold )
    powDosInd <- mean( sapply(res,function(x) x$Dosind$coefficients[11] ) < par$threshold )
    
    powF <- mean( sapply(res,function(x) x$PPf$test[3] ) < par$threshold )
    powDosF <- mean( sapply(res,function(x) x$Dosf$coefficients[11] ) < par$threshold )
    
    powG <- mean( sapply(res,function(x) x$geno$coefficients[11] ) < par$threshold )
    
    diffF<-mean( sapply(res,function(x) exp(x$PPf$full$par[1])-exp(beta) ) )
    diffInd<- mean( sapply(res,function(x) exp(x$PPind$full$par[1])-exp(beta) ) )
    
    diffDosF<-mean( sapply(res,function(x) exp(x$Dosf$coefficients[2])-exp(beta) ) )
    diffDosInd<- mean( sapply(res,function(x) exp(x$Dosind$coefficients[2])-exp(beta) ) )
    
    diffGeno<-mean( sapply(res,function(x) exp(x$geno$coefficients[2])-exp(beta) ) )
    
    betaF<-mean(sapply(res,function(x) x$PPf$full$par[1]))
    betaDosF<-mean(sapply(res,function(x) x$Dosf$coefficients[2]))
    betaGeno<-mean(sapply(res,function(x) x$geno$coefficients[2]))
    
    alphaF<-mean(sapply(res,function(x) x$PPf$full$par[2]))
    alphaDosF<-mean(sapply(res,function(x) x$Dosf$coefficients[1]))
    alphaGeno<-mean(sapply(res,function(x) x$geno$coefficients[1]))
    
    f<- sapply(res,function(x) x$f )
    
  }
  
  return(list(powInd=powInd,powDosInd=powDosInd,powF=powF,powDosF=powDosF,powG=powG,
              diffF=diffF,diffInd=diffInd,diffGeno=diffGeno,diffDosF=diffDosF,diffDosInd=diffDosInd,
              betaF=betaF,betaDosF=betaDosF,betaGeno=betaGeno,alphaF=alphaF,alphaDosF=alphaDosF,alphaGeno=alphaGeno,f=f))
  
}

cex<-1.5

###########################################################

param <- list(sd=1,Nind=1000,beta=0.6,freq=c(0.9,0.1),delta=3,simulations=Nsimulations)
param$Q<-rbind(rep(c(0.05,0.95,0.5,0.95),each=param$Nind/4),1-rep(c(0.05,0.95,0.5,0.95),each=param$Nind/4))
param$fInd <- colSums(param$Q*param$freq)
param$depthGroup<- c(4,1)
param$threshold <- 1e-3

deltaValue<-5
gammaValue<-1


###########################################################

range <- seq(0,deltaValue,length.out=11)

power <- lapply(range,function(x) fun2(delta=x,beta=0,gamma=gammaValue,keepInformative = T,estimateMaf = 0,par=param))
power2 <- lapply(range,function(x) fun2(delta=x,beta=0,gamma=gammaValue,keepInformative = F,estimateMaf = 0,par=param))
power3 <- lapply(range,function(x) fun2(delta=x,beta=0,gamma=gammaValue,keepInformative = T,estimateMaf = 2,par=param))
power4 <- lapply(range,function(x) fun2(delta=x,beta=0,gamma=gammaValue,keepInformative = F,estimateMaf = 2,par=param))

save(range,power,power2,power3,power4,file = "Fig2.Rdata")


######################################
## real power
#########################################

cex<-1.4

range <- seq(0,0.6,length.out=11)

power <- lapply(range,function(x) fun2(delta=0,beta=x,gamma=gammaValue,keepInformative = T,estimateMaf = 0,par=param))
power2 <- lapply(range,function(x) fun2(delta=0,beta=x,gamma=gammaValue,keepInformative = F,estimateMaf = 0,par=param))
power3 <- lapply(range,function(x) fun2(delta=0,beta=x,gamma=gammaValue,keepInformative = T,estimateMaf = 2,par=param))
power4 <- lapply(range,function(x) fun2(delta=0,beta=x,gamma=gammaValue,keepInformative = F,estimateMaf = 2,par=param))


save(range,power,power2,power3,power4,file = "Fig3.Rdata")



########################################
##### SCENARIOS WITHOUT ADMIXTURE
########################################


param <- list(sd=1,Nind=1000,beta=0.6,freq=c(0.45,0.45),delta=3,simulations=Nsimulations)
param$Q<-rbind(rep(c(1,1,1,1),each=param$Nind/4),(1-rep(c(1,1,1,1),each=param$Nind/4)))
param$fInd <- colSums(param$Q*param$freq)
param$depthGroup<- c(4,1)
param$threshold <- 1e-3

deltaValue<-5
gammaValue<-0


###########################
## effect of depth
#################################


range <- seq(0,deltaValue,length.out=11)

## changed from keepInformative = T to keepInformative = F - 10-04-2020
power <- lapply(range,function(x) fun2(delta=x,beta=0,gamma=gammaValue,keepInformative = T,estimateMaf = 0,par=param))
power2 <- lapply(range,function(x) fun2(delta=x,beta=0,gamma=gammaValue,keepInformative = F,estimateMaf = 0,par=param))
power3 <- lapply(range,function(x) fun2(delta=x,beta=0,gamma=gammaValue,keepInformative = T,estimateMaf = 2,par=param))
power4 <- lapply(range,function(x) fun2(delta=x,beta=0,gamma=gammaValue,keepInformative = F,estimateMaf = 2,par=param))



save(range,power,power2,power3,power4,file = "SuppFig1.Rdata")

###############################



range <- seq(0,0.6,length.out=11)

power <- lapply(range,function(x) fun2(delta=0,gamma=gammaValue,beta=x,keepInformative = T,estimateMaf = 0,par=param))
power2 <- lapply(range,function(x) fun2(delta=0,gamma=gammaValue,beta=x,keepInformative = F,estimateMaf = 0,par=param))
power3 <- lapply(range,function(x) fun2(delta=0,gamma=gammaValue,beta=x,keepInformative = T,estimateMaf = 2,par=param))
power4 <- lapply(range,function(x) fun2(delta=0,gamma=gammaValue,beta=x,keepInformative = F,estimateMaf = 2,par=param))


save(range,power,power2,power3,power4,file = "SuppFig2.Rdata")

########################################
##### SCENARIOS WITHOUT ADMIXTURE
########################################

param <- list(sd=1,Nind=1000,beta=0.6,freq=c(0.45,0.45),delta=3,simulations=Nsimulations)
param$Q<-rbind(rep(c(1,1,1,1),each=param$Nind/4),(1-rep(c(1,1,1,1),each=param$Nind/4)))
param$fInd <- colSums(param$Q*param$freq)
param$depthGroup<- c(4,1)
param$threshold <- 1e-3

deltaValue<-5
gammaValue<-0

range <-  seq(0,0.5,by=0.06)


power <- lapply(range,function(x) fun2(delta=deltaValue,gamma=gammaValue,beta=x,keepInformative = T,estimateMaf = 0,par=param))
power2 <- lapply(range,function(x) fun2(delta=deltaValue,gamma=gammaValue,beta=x,keepInformative = F,estimateMaf = 0,par=param))
power3 <- lapply(range,function(x) fun2(delta=deltaValue,gamma=gammaValue,beta=x,keepInformative = T,estimateMaf = 2,par=param))
power4 <- lapply(range,function(x) fun2(delta=deltaValue,gamma=gammaValue,beta=x,keepInformative = F,estimateMaf = 2,par=param))



save(range,power,power2,power3,power4,file = "SuppFig3and4.Rdata")


###############################################
###############################################

param <- list(sd=1,Nind=1000,beta=0.6,freq=c(0.9,0.1),delta=3,simulations=Nsimulations)
param$Q<-rbind(rep(c(0.05,0.95,0.5,0.95),each=param$Nind/4),1-rep(c(0.05,0.95,0.5,0.95),each=param$Nind/4))
param$fInd <- colSums(param$Q*param$freq)
param$depthGroup<- c(4,1)
param$threshold <- 1e-3

deltaValue<-5
gammaValue<-1




range <- seq(0,0.6,length.out=11)

power <- lapply(range,function(x) fun2(delta=deltaValue,beta=x,gamma=gammaValue,keepInformative = T,estimateMaf = 0,par=param))
power2 <- lapply(range,function(x) fun2(delta=deltaValue,beta=x,gamma=gammaValue,keepInformative = F,estimateMaf = 0,par=param))
power3 <- lapply(range,function(x) fun2(delta=deltaValue,beta=x,gamma=gammaValue,keepInformative = T,estimateMaf = 2,par=param))
power4 <- lapply(range,function(x) fun2(delta=deltaValue,beta=x,gamma=gammaValue,keepInformative = F,estimateMaf = 2,par=param))
power5 <- lapply(range,function(x) fun2(delta=deltaValue,beta=x,gamma=gammaValue,keepInformative = T,estimateMaf = 1,par=param))
power6 <- lapply(range,function(x) fun2(delta=deltaValue,beta=x,gamma=gammaValue,keepInformative = F,estimateMaf = 1,par=param))



save(range,power,power2,power3,power4,file = "SuppFig5and6.Rdata")

############################################################
############################################################



param <- list(sd=1,Nind=1000,beta=0.6,freq=c(0.9,0.1),delta=3,simulations=Nsimulations)
param$Q<-rbind(rep(c(0.05,0.95,0.5,0.95),each=param$Nind/4),1-rep(c(0.05,0.95,0.5,0.95),each=param$Nind/4))
param$fInd <- colSums(param$Q*param$freq)
param$depthGroup<- c(4,1)
param$threshold <- 1e-3

deltaValue<-5
gammaValue<-1


range <- seq(0,0.6,length.out=11)

power <- lapply(range,function(x) funBin2(delta=deltaValue,beta=x,gamma=gammaValue,keepInformative = T,estimateMaf = 0,par=param))
power2 <- lapply(range,function(x) funBin2(delta=deltaValue,beta=x,gamma=gammaValue,keepInformative = F,estimateMaf = 0,par=param))
power3 <- lapply(range,function(x) funBin2(delta=deltaValue,beta=x,gamma=gammaValue,keepInformative = T,estimateMaf = 2,par=param))
power4 <- lapply(range,function(x) funBin2(delta=deltaValue,beta=x,gamma=gammaValue,keepInformative = F,estimateMaf = 2,par=param))


save(range,power,power2,power3,power4,file = "SuppFig7and8.Rdata")

###################################################################
##################################################################

RListP<-list()
snptestListP<-list()
snptestListP2<-list()
angsdListP<-list()
angsdListP2<-list()
genoListP<-list()
doseListP<-list()
doseListP2<-list()


betaRange<-seq(0,0.5,length=21)

threshold<-1e-03

for(beta in betaRange){
  
  beta<-beta
  
  Nind<-10000
  
  freq<-c(0.9,0.1)
  ##dep<-c(rpois(Nind/4,2),rpois(Nind/4,2),rpois(Nind/4,1),rpois(Nind/4,1))
  dep <- rep(c(0.1,1,10,20),c(Nind/4,Nind/4,Nind/4,Nind/4))
  q <- rep(c(0.15,0.4,0.95,1),each=Nind/4)
  q <- rbind(q,1-q)
  fInd <- colSums(q*freq)
  
  gamma <- 0.5
  delta<-0.0
  
  print(beta)
  
  for(j in 1:Nsimulations){
    
    freq<-c(0.9,0.1)
    ##dep<-c(rpois(Nind/4,2),rpois(Nind/4,2),rpois(Nind/4,1),rpois(Nind/4,1))
    dep <- rep(c(0.1,1,10,20),c(Nind/4,Nind/4,Nind/4,Nind/4))
    q <- rep(c(0.15,0.4,0.95,1),each=Nind/4)
    q <- rbind(q,1-q)
    fInd <- colSums(q*freq)
    
    gamma <- 0.5
    delta<-0.0
    
    print(j)
    
    g <- rbinom(Nind,2,fInd)
    ##f <- sum(g)/Nind/2
    
    y <- rnorm(Nind)+beta*g + gamma*q[1,] + dep*delta
    
    likes<-getLikes(g,d=dep,norm=T,e=0.01)
    GL <- likes$res       
    
    maf <- squarem( c( 0.1, 0.1 ), emFreqAdmix, GL = t(GL), anc = q[1,] )$par
    fInd <- q[1,]*maf[1] + ( 1-q[1,] )*maf[2]
    f<-estPem(t(GL))
    
    PPf <- GL*c((1-f)^2,2*f*(1-f),f^2)
    PPf <- PPf/rep(colSums(PPf),each=3)
    Dosf<- colSums(PPf*c(0,1,2))
    
    PPind <- GL*rbind((1-fInd)^2,2*fInd*(1-fInd),fInd^2)
    PPind <- PPind/rep(colSums(PPind),each=3)
    Dosind<- colSums(PPind*c(0,1,2))
    
    
    Nkeep<-Nind
    
    ## ANGSD
    write(paste("marker allele1 allele2 ",paste(paste0(rep("Ind",Nkeep*3),rep(1:Nkeep,each=3)),collapse = " ")),"simulated.beagle")
    write(paste("1_1000 A G",paste(as.vector(PPf),collapse = " ")),"simulated.beagle",append=T)
    
    write(paste("marker allele1 allele2 ",paste(paste0(rep("Ind",Nkeep*3),rep(1:Nkeep,each=3)),collapse = " ")),"simulated2.beagle")
    write(paste("1_1000 A G",paste(as.vector(PPind),collapse = " ")),"simulated2.beagle",append=T)
    
    write.table(y,"simulated.phe",col=F,row=F,qu=F)
    write.table(q[1,],"simulated.cov",col=F,row=F,qu=F)
    
    ## SNPTEST
    write(paste("SNP1 1_1000 1000 A G",paste(as.vector(PPf),collapse = " ")),"simulated.gen")
    write(paste("SNP1 1_1000 1000 A G",paste(as.vector(PPind),collapse = " ")),"simulated2.gen")
    
    write(paste("ID_1 ID_2 missing cov_1 pheno1"),"simulated.sample1")
    write(paste("0 0 0 C P"),"simulated.sample1",append = T)
    
    write.table(cbind(1:Nkeep,1:Nkeep,0,q[1,],y),"simulated.sample2",col=F,row=F,qu=F)
    
    system("cat simulated.sample1 simulated.sample2 > simulated.sample")
    
    system("./angsd -doMaf 4 -beagle simulated.beagle -yQuant simulated.phe -doAsso 4 -out simulated -fai hg19.fa.fai -cov simulated.cov")
    system("./angsd -doMaf 4 -beagle simulated2.beagle -yQuant simulated.phe -doAsso 4 -out simulated2 -fai hg19.fa.fai -cov simulated.cov")
    
    system("./snptest_v2.5.4-beta3 -data simulated.gen simulated.sample -frequentist 1 -method em -o simulated.res -pheno pheno1 -cov_all -use_raw_phenotypes -use_raw_covariates ")
    system("./snptest_v2.5.4-beta3 -data simulated2.gen simulated.sample -frequentist 1 -method em -o simulated2.res -pheno pheno1 -cov_all -use_raw_phenotypes -use_raw_covariates ")
    
    PPf2<-EMasso(GP=t(PPf),ys=y,covs=cbind(1,q[1,]),qu=T,logEM = T)
    
    
    
    snptest<-read.table("simulated.res",as.is=T,h=T)
    snptest2<-read.table("simulated2.res",as.is=T,h=T)
    angsd<-read.table("simulated.lrt0.gz",as.is=T,h=T)
    angsd2<-read.table("simulated2.lrt0.gz",as.is=T,h=T)
    
    RListP[[paste0(beta,"_",j)]]<-PPf2$test[3]
    snptestListP[[paste0(beta,"_",j)]]<-snptest$frequentist_add_pvalue
    angsdListP[[paste0(beta,"_",j)]]<-angsd$LRT
    snptestListP2[[paste0(beta,"_",j)]]<-snptest2$frequentist_add_pvalue
    angsdListP2[[paste0(beta,"_",j)]]<-angsd2$LRT
    
    
    RList[[paste0(beta,"_",j)]]<-PPf2$full$par[1]
    snptestList[[paste0(beta,"_",j)]]<-snptest$frequentist_add_beta_1
    angsdList[[paste0(beta,"_",j)]]<-angsd$beta
    snptestList2[[paste0(beta,"_",j)]]<-snptest2$frequentist_add_beta_1
    angsdList2[[paste0(beta,"_",j)]]<-angsd2$beta
    
    
    s<-summary(lm(y~g+q[1,]))
    genoList[[paste0(beta,"_",j)]]<-s$coefficients[2,1]
    genoListP[[paste0(beta,"_",j)]]<-s$coefficients[2,4]
    
    d<-summary(lm(y~Dosf+q[1,]))
    doseList[[paste0(beta,"_",j)]]<-d$coefficients[2,1]
    doseListP[[paste0(beta,"_",j)]]<-d$coefficients[2,4]
    
    d2<-summary(lm(y~Dosind+q[1,]))
    doseList2[[paste0(beta,"_",j)]]<-d2$coefficients[2,1]
    doseListP2[[paste0(beta,"_",j)]]<-d2$coefficients[2,4]
    
    
  }
  
}



genoVector<-sapply(betaRange, function(x) mean(unlist(genoList[grepl(paste0("^",x,"_"),names(genoList))])-as.numeric(x)))

snptestVector<-sapply(betaRange, function(x) mean(unlist(snptestList[grepl(paste0("^",x,"_"),names(snptestList))])-as.numeric(x)))
snptestVector2<-sapply(betaRange, function(x) mean(unlist(snptestList2[grepl(paste0("^",x,"_"),names(snptestList2))])-as.numeric(x)))

angsdVector<-sapply(betaRange, function(x) mean(unlist(angsdList[grepl(paste0("^",x,"_"),names(angsdList))])-as.numeric(x)))
angsdVector2<-sapply(betaRange, function(x) mean(unlist(angsdList2[grepl(paste0("^",x,"_"),names(angsdList2))])-as.numeric(x)))

doseVector<-sapply(betaRange, function(x) mean(unlist(doseList[grepl(paste0("^",x,"_"),names(doseList))])-as.numeric(x)))
doseVector2<-sapply(betaRange, function(x) mean(unlist(doseList2[grepl(paste0("^",x,"_"),names(doseList2))])-as.numeric(x)))


genoVectorP<-sapply(betaRange, function(x) mean(unlist(genoListP[grepl(paste0("^",x,"_"),names(genoListP))])<threshold))

snptestVectorP<-sapply(betaRange, function(x) mean(unlist(snptestListP[grepl(paste0("^",x,"_"),names(snptestListP))])<threshold))
snptestVectorP2<-sapply(betaRange, function(x) mean(unlist(snptestListP2[grepl(paste0("^",x,"_"),names(snptestListP2))])<threshold))

angsdVectorP<-sapply(betaRange, function(x) mean(unlist(angsdListP[grepl(paste0("^",x,"_"),names(angsdListP))])>qchisq(1-threshold,df=1)))
angsdVectorP2<-sapply(betaRange, function(x) mean(unlist(angsdListP2[grepl(paste0("^",x,"_"),names(angsdListP2))])>qchisq(1-threshold,df=1)))

doseVectorP<-sapply(betaRange, function(x) mean(unlist(doseListP[grepl(paste0("^",x,"_"),names(doseListP))])<threshold))
doseVectorP2<-sapply(betaRange, function(x) mean(unlist(doseListP2[grepl(paste0("^",x,"_"),names(doseListP2))])<threshold))


save(file="SuppFig12.Rdata",betaRange,genoVector,snptestVector,snptestVector2,doseVector,doseVector2,
     angsdVector,angsdVector2,genoVectorP,snptestVectorP,snptestVectorP2,angsdVectorP,angsdVectorP2,doseVectorP,doseVectorP2)
