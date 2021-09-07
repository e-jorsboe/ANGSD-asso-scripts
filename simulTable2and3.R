
library(parallel)

source("ANGSD-asso/regression_from_design.R")
source("ANGSD-asso/ANGSD-assoR.R")
source("ANGSD-asso/ANGSD-assoSimR.R")


args<-commandArgs(trailingOnly = T)

cases<-as.numeric(args[1])
controls<-as.numeric(args[2])

depCases<-as.numeric(args[3])
depControls<-as.numeric(args[4])

RR<-as.numeric(args[5])

keepAllanalysis<-as.numeric(args[6])

cores<-as.numeric(args[7])

effect <- seq(1,2,by=0.1)
dep<-c(rep(depCases,cases),rep(depControls,controls))

simulations<-10000

freq<-0.05
yBin<-c(rep(1,cases),rep(0,controls))

binScore<-function(y,GP){
  if(is.null(dim(GP))){
    mat<-matrix(0,nrow=3,ncol=length(y))
    mat[GP+(1:length(y)-1)*3+1]<-1
    GP<-t(mat)
  }
  
  Ex<-as.vector(GP%*%0:2)
  Ex2<-as.vector(GP%*%(0:2)^2)
  yTilde<-mean(y)
  
  U<-sum((y-yTilde)*Ex)
  Vaa<-yTilde*(1-yTilde)*length(y)
  Vab<-yTilde*(1-yTilde)*sum(Ex)
  Vbb<-sum((yTilde*(1-yTilde)-(y-yTilde)^2)*Ex2+(y-yTilde)^2*Ex^2)
  I<-Vbb-Vab^2/Vaa
  1-pchisq(U^2/I,1)
}


normScore<-function(y,GP){
  if(is.null(dim(GP))){
    mat<-matrix(0,nrow=3,ncol=length(y))
    mat[GP+(1:length(y)-1)*3+1]<-1
    GP<-t(mat)
  }
  
  Ex<-as.vector(GP%*%0:2)
  Ex2<-as.vector(GP%*%(0:2)^2)
  yTilde<-mean(y)
  var<-var(y)
  
  U<-sum((y-yTilde)/var*Ex)
  Vaa<-1/var*length(y)
  Vab<-1/var*sum(Ex)
  k<-(y-yTilde)^2/var^2
  Vbb<-sum( (1/var-k)*Ex2+k*Ex^2)
  I<-Vbb-Vab^2/Vaa
  1-pchisq(U^2/I,1)
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


## rr is relative risk, p is maf, kp is prevalence of disease 
freq.rr<-function (rr = 4.5, p = 0.15, kp = 0.1,freq=F){#kp prevalence
  f.mod <- c(rr^2, rr,1)
  q <- 1 - p
  fhw <- c(p^2, 2 * p * q, q^2)
  pi <- kp/sum(f.mod * fhw)
  if (pi <= 0 | pi >= 1) {
    warning("The combination of p, kp, and rr produces an unrealistic value of pi.")
    ret <- NA
  }
  else {
    fe <- rbind(fhw, fhw)
    dimnames(fe) <- list(c("Case", "Control"), c("AA", "Aa","aa"))
    f <- fe * rbind(f.mod * pi, 1 - f.mod * pi)
    ret<-f/ apply(f, 1, sum)
    if(freq)
      ret<-(ret[,2]+2*ret[,1])/2
  }
  ret
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


funBin <- function(seed,beta,maf,prevalence,keepAll=FALSE){
  set.seed(seed)  
  freqs<-freq.rr(rr=beta,p=maf,kp=prevalence,freq = T)  
  g <- c(rbinom(cases,2,freqs[1]),rbinom(controls,2,freqs[2]))
  likes<-getLikes(g,d=dep,norm=T,e=0.01)
  GL <- likes$res  
  keep<-likes$depth>0
  
  if(keepAll){
    ## estimate freq using EM algo
    f<- estPem(t(GL))
    
    PPf <- GL*c((1-f)^2,2*f*(1-f),f^2)
    
    PPf <- PPf/rep(colSums(PPf),each=3)
    Dosf<- colSums(PPf*c(0,1,2))
    
    g2<-g
    yBin2<-yBin
    
  } else{
  
    ## estimate freq using EM algo
    f<- estPem(t(GL[,keep]))
    
    PPf <- GL[,keep]*c((1-f)^2,2*f*(1-f),f^2)
    
    PPf <- PPf/rep(colSums(PPf),each=3)
    Dosf<- colSums(PPf*c(0,1,2))
    
    g2<-g[keep]
    yBin2<-yBin[keep]
    
  }
    
    fcases<-sum(Dosf[yBin2==1]**2)/(2*length(Dosf[yBin2==1]))
    fcontrols<-sum(Dosf[yBin2==0]**2)/(2*length(Dosf[yBin2==0]))
    
    r2cases<-(sum(Dosf[yBin2==1]**2)/length(Dosf[yBin2==1])-(sum(Dosf[yBin2==1])/length(Dosf[yBin2==1]))**2) / (2*fcases*(1-fcases))
    r2controls<-(sum(Dosf[yBin2==0]**2)/length(Dosf[yBin2==0])-(sum(Dosf[yBin2==0])/length(Dosf[yBin2==0]))**2) / (2*fcontrols*(1-fcontrols))
    
    s<-apply(PPf,2,function(x) max(x)-min(x))
    
    list(
        score=binScore(y=yBin2,GP=t(PPf)),
        PPf=ANGSDasso(GP=t(PPf),ys=yBin2,qu=F),
        Dosf=summary(glm(yBin2~Dosf,family="binomial")),
        geno=summary(glm(yBin2~g2,family="binomial")),
        r2cases=r2cases,
        r2controls=r2controls
    )
}


funBin2 <- function(beta,maf,prevalence,keepAll=FALSE){
  cat(" beta",beta,"\n")
  res <- parallel::mclapply(1:simulations, function(x) funBin(x,beta=beta,maf=maf,prevalence = prevalence,keepAll=keepAll ),mc.cores=cores)
  
  powF <- mean( sapply(res,function(x) x$PPf$test[3] ) < 0.00001 )
  powDosF <- mean( sapply(res,function(x) x$Dosf$coefficients[8] ) < 0.00001 )
 
  powG <- mean( sapply(res,function(x) x$geno$coefficients[8] ) < 0.00001 )
  
  score <- mean( sapply(res,function(x) x$score ) < 0.00001 )
  
  diffF<-mean( sapply(res,function(x) x$PPf$full$par[1]-x$geno$coefficients[2] ) )
  
  diffDosF<-mean( sapply(res,function(x) x$Dosf$coefficients[2]-x$geno$coefficients[2] ) )

  r2cases <- mean( sapply(res,function(x) x$r2cases ))
  r2controls <- mean( sapply(res,function(x) x$r2controls ))
  
  c(score=score,powerF=powF,powDosF=powDosF,powerG=powG,diffF=diffF,diffDosF=diffDosF,r2cases=r2cases,r2controls=r2controls)
}



if(keepAllanalysis>0){
    
    a0<-funBin2(beta=RR,maf=0.05,prevalence = 0.1,keepAll=TRUE)
    
    save(a0,file=paste0("casesControl_RR",RR,"_Nind",(controls+cases),"depth",depCases,"_",depControls,".Rdata"))
                    
    
} else{
    
    f0<-funBin2(beta=RR,maf=0.05,prevalence = 0.1)
    
    save(f0,file=paste0("casesControl_RR",RR,"_Nind",(controls+cases),"depth",depCases,"_",depControls,".rmDepth0.Rdata"))        
    
}
