## x=genotype
##d=dybde (avg)
## e=error rate
## norm=Normalize to sum to one
getLikes<-function(x,d=2,e=0.01,norm=FALSE){
    n<-length(x)
    dep<-rpois(n,d)
    nA<-rbinom(n,dep,c(e,0.5,1-e)[x+1])
    res<-rbind(dbinom(nA,dep,e),dbinom(nA,dep,0.5),dbinom(nA,dep,1-e))
    if(norm){
        res<-t(t(res)/colSums(res))
    }
    res[is.na(res)]<-1/3
    return(res)
}

## function for converting beagle file to matrix with GPs for EMasso
beagleTo<-function(beagle){
    cols<-ncol(beagle)
    rows<-ncol(beagle)    
    m<-matrix(unlist(beagle[,4:cols]),ncol=3,nrow=((rows-3)/3),byrow = T)
    return(m)    
}



## Simulate ancestry specific assocation
##
## param qvec individual admixture proportions from population 1 (n vec)
## param pvec allele frequencies in the two ancestral populations (2 vec)
## param bvec population specific effects (2 vec for add and 3 vec for rec) or (4 long if you want to simulate effect for un-counted allele)
## param model either "add" or "rec"
## param cvec intercept and effects of covariates if any (p vec)
## param cov design matrix with intercept (p vec)
## param ancE effect of ancestry in SD, by default no effect (0,0)
## param quant simulate quantitative trait, logistic regression model if FALSE
##return list of genotype vector (gt[n]), phenotype vector (y[n]), ancestral state (as[2,n]), allelic genotypes (ag[2,n]) and design for simulated population specific effects (x [2,n] for additive and [3,n] for recessive)
ANGSDassoSim <- function(qvec, pvec, bvec, model, cvec, cov, ancE=c(0,0),quant=TRUE,setSeed=NULL,depth,depthEffect=0){

    if(!is.null(setSeed)){
        set.seed(setSeed)
    }

    nInd <- length(qvec)
    ## add intersect
    if(missing(cov)){
        cov <- matrix(1, nrow=nInd, ncol=1)
    }

    qmat <- rbind(qvec,1-qvec)
    ## Allelic state (1=pop1, 2=pop2)
    as<-apply(qmat, 2, function(x)sample(1:2,2,replace=TRUE,prob=x)) ## 2 x nInd
    ## Sample allelic genotype (1=risk allele)
    ag <- matrix(rbinom(2*nInd,1,pvec[as]),nrow=2,byrow=FALSE) ## 2 x nInd

    ## checking that bvec is either 2 or 4 long for the additive model
    if(model=="add" & length(bvec)!=2){
        stop("bvec must have length of 2")
    }

    ## checking that bvec is either 3 or 5 long for the recessive model
    if(model=="rec" & length(bvec)!=3){
        stop("bvec must have length of 3")
    }

    ## simulating states for the additive model
    if(model=="add") {
        x1 <- as.integer(ag[1,]*(as[1,]==1)+ag[2,]*(as[2,]==1))
        x2 <- as.integer(ag[1,]*(as[1,]==2)+ag[2,]*(as[2,]==2))
        x <- cbind(x1,x2)
    } else if(model=="rec") {
        x1 <- as.integer((colSums(as) == 2) * (colSums(ag) == 2))
        x2 <- as.integer((colSums(as) == 4) * (colSums(ag) == 2))
        xM <- as.integer((colSums(as) == 3) * (colSums(ag) == 2))
        x <- cbind(x1,x2,xM)
    } else {
        stop('Model must be either "add" or "rec", please.')
    }

    if(model=="add"){
        ## also incorporating effects from both alternative allele (pop 1 and 2)
        evec <- x %*% bvec + as.matrix(cov) %*% cvec + ifelse(as[1,]==1,ancE[1],ancE[2]) + ifelse(as[2,]==1,ancE[1],ancE[2]) + depth*depthEffect
    } else{
        ## also incorporating effects from both recessive alternative allele (pop 1 and 2)
        evec <- x %*% bvec + as.matrix(cov) %*% cvec + ifelse(as[1,]==1,ancE[1],ancE[2]) + ifelse(as[2,]==1,ancE[1],ancE[2]) + depth*depthEffect
    }
    
    ## Simulated responce
    ## note that standard deviation is 1
    if(quant){
        yvec<-rnorm(n = nInd, mean = evec)
    }else{
        yvec <- rbinom(n = nInd, size = 1, prob = plogis(evec))
    }
    gt<- colSums(ag)

    GL<-as.vector(getLikes(gt,d=depth,norm=T))

    return(list(gt = gt, y = yvec, as = as, ag = ag, x = x,evec=evec,GL=GL))
}


gl<-ANGSDassoSim(qvec =  rep(c(0,0.25,0.5,0.75,1),each=2500/5) ,pvec = c(0.2,0.2),bvec=c(0,0),ancE=c(0.2,0),model="add",quant = F,cov=rep(1,2500),cvec=0,depth=rep(c(1,3,5,10,20),each=2500/5))


tmp<-ANGSDassoSim(qvec =  rep(c(0,0.25,0.5,0.75,1),each=2500/5) ,pvec = c(0.1,0.3), bvec=c(0,0),ancE=c(0,0.3),model="add",quant = F,cov=rep(1,2500),cvec=0,depth=rep(c(1,2,3,4,5),each=2500/5))
tmp2<-ANGSDassoSim(qvec =  rep(c(0,0.25,0.5,0.75,1),each=2500/5) ,pvec = c(0.1,0.3), bvec=c(0,0),ancE=c(0,0.3),model="add",quant = F,cov=rep(1,2500),cvec=0,depth=rep(c(1,2,3,4,5)/10,each=2500/5))
