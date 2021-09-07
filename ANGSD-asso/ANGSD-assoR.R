


#' Making prior distribution across states
#'
#' @param geno genotype vector (length n)
#' @param admix admixture proportion from pop 1 (vector length n)
#' @param freq allele frequencies in pop 1 and 2 (vector length 2)
#' @return Probability of ancestral state conditional on genotype (vector length 4n with one entry per possible ancestral configuration \{1,1\}, \{1,2\}, \{2,1\}, \{2,2\} for each individual)
geno_prob <- function(GL, freq) {
    GL2<-t(t(GL)*c(freq**2,2*freq*(1-freq),(1-freq)**2))
    GL3<-rowSums(GL2)
    return(GL/GL3)    
}


#' Calculate joint probabilities of phenotypes and ancestral states for quantitative traits
#'
#' @param params Parameters, first genetic effects, then intercept and additional covariate effects and last sd (length 2/3 + p + 1)
#' @param phenos Individual phenotypes (augmented to length 4n, i.e. each repeated 4 times)
#' @param design Augmented design matrix ( matrix 4n x length(params)-1 )
#' @param prior Probability of ancestral state conditional on genotype (vector length 4n) - see ?make_mixture_props
#' @return Joint probability of phenotype and state given parameters (matrix 4 x n)
make_joint_lm <- function(params, phenos, design, prior) {
    
    probs<-dnorm(x = phenos,mean = design%*%(params[-length(params)]),sd = params[length(params)])    
    matrix(probs * prior,nr=3)
}

#' Calculate log of joint probabilities of phenotypes and ancestral states for quantitative traits
#'
#' @param params Parameters (first genetic effects, then intercept and additional covariate effects and last sd)
#' @param phenos Individual phenotypes augmented to length 4 x N
#' @param design Augmented design matrix, 4N rows
#' @return Matrix with 4 rows and N columns, with log of joint probability of phenotype and state given parameters
make_joint_lm_log <- function(params, phenos, design, prior) {
    
    probs <- dnorm(x = phenos, mean = design%*%(params[-length(params)]),
                   sd = params[length(params)], log=T)   
    matrix(probs + log(prior), nr=3)
}


#' Calculate joint probabilities of phenotypes and ancestral states for dichotomuous traits
#'
#' @param params Parameters (first genetic effects, then intercept and additional covariate effects and last sd)
#' @param phenos Individual phenotypes augmented to length 4 x N
#' @param design Augmented design matrix, 4N rows
#' @return Matrix with 4 rows and N columns, with joint probability of phenotype and state given parameters
make_joint_bin <- function(params, phenos, design, prior) {
    
    probs <- dbinom(x=phenos,size=1,prob=plogis(design%*%params))               
    matrix(probs * prior,nr=3)
}

#' Calculate joint probabilities of phenotypes and ancestral states for dichotomuous traits
#'
#' @param params Parameters (first genetic effects, then intercept and additional covariate effects and last sd)
#' @param phenos Individual phenotypes augmented to length 4 x N
#' @param design Augmented design matrix, 4N rows
#' @return Matrix with 4 rows and N columns, with joint probability of phenotype and state given parameters
make_joint_bin_log <- function(params, phenos, design, prior) {

    probs <- dbinom(x=phenos,size=1,prob=plogis(design%*%params),log=T)

    minval<-min(probs[is.finite(probs)])
    maxval<-max(probs[!(probs%in%c(0))])
        
    probs[probs==0]<-maxval
    probs[is.infinite(probs)]<-minval  
        
    matrix(probs + log(prior),nr=3)          
}



#' Calculate joint probabilities of phenotypes and ancestral states for count traits
#'
#' @param params Parameters (first genetic effects, then intercept and additional covariate effects and last sd)
#' @param phenos Individual phenotypes augmented to length 4 x N
#' @param design Augmented design matrix, 4N rows
#' @return Matrix with 4 rows and N columns, with joint probability of phenotype and state given parameters
make_joint_pois <- function(params, phenos, design, prior) {
    
    probs <- dpois(x=phenos,lambda=exp(design%*%params))         
    matrix(probs * prior,nr=3)
}

#' Calculate joint probabilities of phenotypes and ancestral states for count traits
#'
#' @param params Parameters (first genetic effects, then intercept and additional covariate effects and last sd)
#' @param phenos Individual phenotypes augmented to length 4 x N
#' @param design Augmented design matrix, 4N rows
#' @return Matrix with 4 rows and N columns, with joint probability of phenotype and state given parameters
make_joint_pois_log <- function(params, phenos, design, prior) {
    
    probs <- dpois(x=phenos,lambda=exp(design%*%params),log=T)        
    matrix(probs + log(prior),nr=3)
}

#' EM update step for quantitative traits
#'
#' @return list of two: pars_new are new parameter estimates (effects and sd) and mll_new is the updated minus log likelihood
update_em_lm <- function(pars, phen, dsgn, prior,  W, SNPTEST=FALSE) {
  ## joint pheno and state posterior prob
  joint_pheno_state <- make_joint_lm(params = pars, phenos = phen, design = dsgn, prior = prior)
  ## Marginal probs of y (state summed out)
  prob_pheno <- colSums(joint_pheno_state)
  ## minus log likelihood (for old parameters actually)
  mll_new <- -sum(log(prob_pheno))
  ## degrees of freedom for calculatind new sd estimate
  df <- length(prob_pheno) - ncol(dsgn)
  ## Weights are the posterior state probabilities given pheno and pars
  wgts <- c(t(t(joint_pheno_state) / prob_pheno))
  wgts<-wgts*rep(W,each=3)
      
  ## New parameter guess  
  fit <- fit_lm(design = dsgn, pheno = phen, weights = wgts)
  #fit <- lm(phen~.-1,weights = wgts,data = data.frame(phen,dsgn))
  sigm <- sqrt(sum((fit$residuals^2) * wgts) / df)
  pars_new <- c(fit$coefficients, sigm)
  
  if(SNPTEST){
    return(list(pars_new = pars_new, mll_new = mll_new, wgts = wgts))
  } else{
    return(list(pars_new = pars_new, mll_new = mll_new))
  }
}

#' colSums in log-space
colSumsLog <- function(x){
    ## each column in x is log(val1), log(val2), log(val3), log(val4)
    ## so max log value gives me one, which I call max
    k <- apply(x, 2, max)
    ## this is log(exp(log(val1)-log(max)) + exp(log(val2)-log(max)) + exp(log(val3)-log(max))  + exp(log(val4)-log(max))) + log(max)
    ## same (exp(log(val1))/exp(log(max)))*(max)
    log(rowSums(exp(t(x) - k))) + k
}

#' EM update step for quantitative traits, calculations done in log space
#'
#' @param pars Model parameters (length 2/3 + p +1 )
#' @param phen Individual phenotypes augmented to length 4 x N
#' @param dsgn Augmented design matrix (matrix 4n x length(params)-1)
#' @param prior Probability of ancestral state conditional on genotype (vector length 4n) - see ?make_mixture_props
#' @return list of two: pars_new are new parameter estimates (effects and sd) and mll_new is the updated minus log likelihood
update_em_lm_log<-function(pars, phen, dsgn, prior, W) {
    ## joint pheno and state posterior prob
    joint_pheno_state_log <- make_joint_lm_log(params = pars, phenos = phen,
                                               design = dsgn, prior = prior)
    
    prob_pheno_log <- colSumsLog(joint_pheno_state_log)
    mll_new<- -sum(prob_pheno_log)
    df <- length(phen)/3-ncol(dsgn)
    wgts_log<- c(t(t(joint_pheno_state_log) - prob_pheno_log))
    wgts <- exp(wgts_log)    
    wgts<-wgts*rep(W,each=3)
        
    ## New parameter guess
    fit <- fit_lm(design = dsgn, pheno = phen, weights = wgts)
    ##fit <- lm(phen~.-1,weights = wgts,data = data.frame(phen,dsgn))
    sigm <- sqrt(sum((fit$residuals^2) * wgts) / df)
    pars_new <- c(fit$coefficients, sigm)
    
    return(list(pars_new = pars_new, mll_new = mll_new))
}



#' EM update step for dichotomuous traits
#'
#' @return list of two: pars_new are new parameter estimates (effects) and mll_new is the updated minus log likelihood
update_em_bin <- function(pars, phen, dsgn, prior, W) {
  ## joint pheno and state posterior prob
  joint_pheno_state <- make_joint_bin(params = pars, phenos = phen, design = dsgn, prior = prior)
  ## Marginal probs of y (state summed out)
  prob_pheno <- colSums(joint_pheno_state)
  ## minus log likelihood (for old parameters actually)
  mll_new <- -sum(log(prob_pheno))
  ## Weights are the posterior state probabilities given pheno and pars
  wgts <- c(t(t(joint_pheno_state) / prob_pheno))
  wgts<-wgts*rep(W,each=3)
    
  ## New parameter guess
  fit <- fit_bin(design = dsgn, pheno = phen, weights = wgts)
  pars_new <- fit$coefficients
  return(list(pars_new = pars_new, mll_new = mll_new))
}

#' EM update step for dichotomuous traits
#'
#' @return list of two: pars_new are new parameter estimates (effects) and mll_new is the updated minus log likelihood
update_em_bin_log <- function(pars, phen, dsgn, prior, W) {
    ## joint pheno and state posterior prob
    joint_pheno_state_log <- make_joint_bin_log(params = pars, phenos = phen,
                                                 design = dsgn, prior = prior)    
    prob_pheno_log <- colSumsLog(joint_pheno_state_log)    
    mll_new<- -sum(prob_pheno_log)
    
    df <- length(phen)/3-ncol(dsgn)
    wgts_log<- c(t(t(joint_pheno_state_log) - prob_pheno_log))
    wgts <- exp(wgts_log)
    wgts<-wgts*rep(W,each=3)

    ## New parameter guess
    fit <- fit_bin(design = dsgn, pheno = phen, weights = wgts)
  
    pars_new <- c(fit$coefficients)
    return(list(pars_new = pars_new, mll_new = mll_new))
    
}


#' EM update step for dichotomuous traits
#'
#' @return list of two: pars_new are new parameter estimates (effects) and mll_new is the updated minus log likelihood
update_em_pois <- function(pars, phen, dsgn, prior, W) {
  ## joint pheno and state posterior prob
  joint_pheno_state <- make_joint_pois(params = pars, phenos = phen, design = dsgn, prior = prior)
  ## Marginal probs of y (state summed out)
  prob_pheno <- colSums(joint_pheno_state)
  ## minus log likelihood (for old parameters actually)
  mll_new <- -sum(log(prob_pheno))
  ## Weights are the posterior state probabilities given pheno and pars
  wgts <- c(t(t(joint_pheno_state) / prob_pheno))
  wgts<-wgts*rep(W,each=3)
  ## New parameter guess
  fit <- fit_pois(design = dsgn, pheno = phen, weights = wgts)
  pars_new <- fit$coefficients
  return(list(pars_new = pars_new, mll_new = mll_new))
}


#' EM update step for dichotomuous traits
#'
#' @return list of two: pars_new are new parameter estimates (effects) and mll_new is the updated minus log likelihood
update_em_pois_log <- function(pars, phen, dsgn, prior, W) {
    ## joint pheno and state posterior prob
    joint_pheno_state_log <- make_joint_pois_log(params = pars, phenos = phen,
                                                 design = dsgn, prior = prior)
    
    prob_pheno_log <- colSumsLog(joint_pheno_state_log) 
    mll_new<- -sum(prob_pheno_log)

    df <- length(phen)/3-ncol(dsgn)
    wgts_log<- c(t(t(joint_pheno_state_log) - prob_pheno_log))
    wgts <- exp(wgts_log)    
    wgts<-wgts*rep(W,each=3)
   
    ## New parameter guess
    fit <- fit_pois(design = dsgn, pheno = phen, weights = wgts)
  
    pars_new <- c(fit$coefficients)
    return(list(pars_new = pars_new, mll_new = mll_new))
    
}


#' EM update step wrapper for choosing quantitative or dichotomuous traits
update_em_log <- function(pars, phen, dsgn, prior, quant, count, W) {
  if(count == TRUE) {
    #update_em_lm(pars, phen, dsgn, prior)
    update_em_pois_log(pars, phen, dsgn, prior, W)        
  } else if (quant == TRUE) {
    #update_em_lm(pars, phen, dsgn, prior)
    update_em_lm_log(pars, phen, dsgn, prior, W)
  } else {
    update_em_bin_log(pars, phen, dsgn, prior, W)
  }
}

#' EM update step wrapper for choosing quantitative or dichotomuous traits
update_em <- function(pars, phen, dsgn, prior, quant, count, W) {
  if(count == TRUE) {
    #update_em_lm(pars, phen, dsgn, prior)
    update_em_pois(pars, phen, dsgn, prior, W)        
  } else if (quant == TRUE) {
    #update_em_lm(pars, phen, dsgn, prior)
    update_em_lm(pars, phen, dsgn, prior, W)
  } else {
    update_em_bin(pars, phen, dsgn, prior, W)
  }
}

#' Minus log likelihood for full model for quantitative traits (only used by optim checks)
#'
#' @param params Numerical vector of coefficients (and sigma for quantitative)
#' @param desi Augmented design matrix (4 entries per observation) (4n x p)
#' @param phen Augmented phenotypes length 4n
#' @param prior Numerical vector 4 entries per obs of conditional prob of state given genotypes
#' @param quant Logical: is trait quantitiative?
#' @return Minus log of likelihood
calc_mll_lm <- function(params, desi, phen, prior){
  joint_phen_state<-make_joint_lm(params = params, phenos = phen, design = desi, prior = prior)
  return(-sum(log(colSums(joint_phen_state))))
}

#' Minus log likelihood for full model for dicho traits (only used by optim checks)
#'
#' @param params Numerical vector of coefficients (and sigma for quantitative)
#' @param desi Augmented design matrix (4 entries per observation) (4n x p)
#' @param phen Augmented phenotypes length 4n
#' @param prior Numerical vector 4 entries per obs of conditional prob of state given genotypes
#' @param quant Logical: is trait quantitiative?
#' @return Minus log of likelihood
calc_mll_bin <- function(params, desi, phen, prior){
  joint_phen_state <- make_joint_bin(params = params, phenos = phen, design = desi, prior = prior)      
  return(-sum(log(colSums(joint_phen_state))))
}


#' Minus log likelihood for full model for count traits (only used by optim checks)
#'
#' @param params Numerical vector of coefficients (and sigma for quantitative)
#' @param desi Augmented design matrix (4 entries per observation) (4n x p)
#' @param phen Augmented phenotypes length 4n
#' @param prior Numerical vector 4 entries per obs of conditional prob of state given genotypes
#' @param quant Logical: is trait quantitiative?
#' @return Minus log of likelihood
calc_mll_pois <- function(params, desi, phen, prior){
  joint_phen_state <- make_joint_pois(params = params, phenos = phen, design = desi, prior = prior)      
  return(-sum(log(colSums(joint_phen_state))))
}



#' EM algorithm controller
#'
#' @param initial Start values for optimization (if quantitative length is ncol(desi)+1 else length is ncol(desi))
#' @param maxI Max number iterations of EM algo
#' @param phe Observed phenotypes 4n long
#' @param desi Sesign matrix 4n times 2+nCov
#' @param pri Prior dist over states, 4n long
#' @param qua Is trait quantitative? true/false
#' @param tole Convergence tolerence
#' @return list of par (estimates), value (minus log likelihood), counts (iterations), convergence and about (how did algo terminate)
control_em <- function(initial, maxI, phe, desi, pri, qua, tole, count, logEM, W){
    pars_old <- initial
    mll_old <- Inf
    for (iter in 1:maxI) {
         if(logEM){
            update <- update_em_log(pars = pars_old, phen = phe, dsgn = desi, prior = pri, quant=qua, count=count, W=W)
        } else{
            update <- update_em(pars = pars_old, phen = phe, dsgn = desi, prior = pri, quant=qua, count=count, W=W)
        }
        pars_new <- update$pars_new
        mll_new <- update$mll_new
          
        if (mll_new > mll_old + 10e-6) {
            em_out <- list(par=pars_new, value = mll_new, counts = c(iter, iter), convergence = 1, about = "mll_new > mll_old + 10e-6")
            ##print("EM step in wrong direction.")
            break
        } else if (mll_old - mll_new < tole & !is.na(mll_old - mll_new)) {
            em_out <- list(par = pars_new, value = mll_new, counts = c(iter, iter), convergence = 0, about = "mll_old - mll_new < tol")
            break
        } else if ( iter == maxI){
            em_out <- list(par = pars_new, value = mll_new, counts = c(iter, iter), convergence = 1, about="iter == maxI")
        } else {
            pars_old <- pars_new
            if(is.infinite(mll_new)){                
                mll_old <- 999999
            } else{
                mll_old <- mll_new
            }
        }
    }
    
    return(em_out)
}


#' Calculating score for getting observed information quantitative traits
#'
calc_score_vector_lm <- function(esti, phen, dsgn, prior) {
  ## Joint distribution of phenotypes and states
  ##print(esti)
  joint_pheno_state <- make_joint_lm(params = esti, phenos = phen, design = dsgn, prior = prior)
  ## Marginal post probs of y (state summed out) - used for normalizing below
  prob_pheno <- colSums(joint_pheno_state)
  ## Weights are the posterior state probabilities given pheno and pars
  wgts <- c(t(t(joint_pheno_state) / prob_pheno))
  etas <- c(dsgn %*% (esti[-length(esti)]))
  std_resid_wgt <- ((phen-etas) / (esti[length(esti)])^2) * wgts
  matrix(colSums(matrix(c(dsgn * std_resid_wgt), nrow=3)), nrow=length(prob_pheno))
}


#' Calculating score for getting observed information dichotomuous traits
#'
calc_score_vector_bin <- function(esti, phen, dsgn, prior) {
    ## Joint distribution of phenotypes and states
    joint_pheno_state <- make_joint_bin(params = esti, phenos = phen, design = dsgn, prior = prior)
    ## Marginal post probs of y (state summed out) - used for normalizing below
    prob_pheno <- colSums(joint_pheno_state)
    ## Weights are the posterior state probabilities given pheno and pars
    wgts <- c(t(t(joint_pheno_state) / prob_pheno))
    ## etas <- c(dsgn %*% esti)
    etas <- exp(c(dsgn %*% esti))/(1+exp(c(dsgn %*% esti)))
    std_resid_wgt <- (phen-etas) * wgts
    matrix(colSums(matrix(c(dsgn * std_resid_wgt), nrow=3)), nrow=length(prob_pheno))
}

#' Calculating score for getting observed information dichotomuous traits
#'
calc_score_vector_pois <- function(esti, phen, dsgn, prior) {
  ## Joint distribution of phenotypes and states
  joint_pheno_state <- make_joint_pois(params = esti, phenos = phen, design = dsgn, prior = prior)
  ## Marginal post probs of y (state summed out) - used for normalizing below
  prob_pheno <- colSums(joint_pheno_state)
  ## Weights are the posterior state probabilities given pheno and pars
  wgts <- c(t(t(joint_pheno_state) / prob_pheno))
  ## etas <- c(dsgn %*% esti)
  etas <- exp(c(dsgn %*% esti))
  std_resid_wgt <- (phen-etas) * wgts
  matrix(colSums(matrix(c(dsgn * std_resid_wgt), nrow=3)), nrow=length(prob_pheno))
}



#' Calculating score for getting observed information matrix
  calc_score_vector <- function(esti, phen, dsgn, prior, quan, count) {
  if (count) {
    return(calc_score_vector_pois(esti, phen, dsgn, prior))
  } else if (quan) {
    return(calc_score_vector_lm(esti, phen, dsgn, prior))
  } else {
    return(calc_score_vector_bin(esti, phen, dsgn, prior))
  }
}


#' Fit ancestry specific models and test different hypothesis (done after removing obs with missing data!)
#'
#' @param GP matrix N x 3
#' @param ys Phenotype vector of length N
#' @param qvec Admixture proportions from population 1 vector length N
#' @param mafs Frequency of the counted alleles in pop1 and pop2 vector length N
#' @param ancprobs Probabilities of ancestry for each haplotype, matrix (3 cols x N rows) OR (4 cols x N rows) 
#' @param covs Covariate matrix. If missing then intercept is modelled else if intercept is desired then its design matrix must be included in covs matrix.
ANGSDasso <- function(GP, ys, # mandatory arguments 
                   covs, start, # can be missing
                   model='add', qu=TRUE, count=FALSE, maxIter=200, tol=10e-5, ## emil changed maxIter to 200
                   se=FALSE, logEM=FALSE, SNPTEST=FALSE, depths=NULL) {
    if(missing(covs))
        covs<-matrix(1,nrow=length(ys),ncol=1)
    
    if(ncol(GP)!=3 | nrow(GP)!=length(ys)){
      print("not genotype probs - has to have 3 cols and N rows")
      stop()
    }

    if(qu & count){
      print("have to specify either count or quantitative data")
      stop()
    }

    ## ASSUMES THAT INTERCEPT IS PART OF COVARIATES
    if(!all(covs[,1]==1)){
      print("seems like there might not be an intercept in covariates - HAS TO BE INCLUDED!!")
    }
    
    W<-rep(1,length(ys))
    
    if(!is.null(depths)){
        
        if(length(depths)==length(ys)){
          W<-depths            
        } else{
          print("have to be depth for each individual - vector has to be as long as pheno")
          stop()          
        }
    }
     
    ## Augmented dataset 3 times size of original
    pheno <- rep(ys,each=3)
    
    n<-length(ys)
    
    if(model=='add'){    
        design<-cbind(rep(c(0,1,2),times=n),covs[rep(1:n,each=3),])
    } else if(model=='rec'){
        design<-cbind(rep(c(0,0,1),times=n),covs[rep(1:n,each=3),])
    } else{
        print("not valid model")
        stop()
    }
    
    if(missing(start)){
        start <- runif(ncol(design),-1,1)     
        if(qu){
            start <- c(start,sd(ys))
        }
    }

    startNull<-start[-1]
    designNull<-as.matrix(design[,-1])
    
    design<-design
    
    if(SNPTEST) {
      ## EMULATING METHOD USED IN SNPTEST's EM METHOD 
        
      ## running null model without genotype - to be able to calculate residuals
      null<-control_em(initial = startNull, phe = pheno, desi = designNull, pri = t(GP),
                       qua = qu, tole = tol, maxI = maxIter, count=count, logEM=logEM, W=W)
                
      ## removing effect of covariates by regressing them out
      residuals<-(ys - covs%*%null$par[1:(length(null$par)-1)])        
      residuals2 <- rep(residuals,each=3)

      ## then running EM model on residuals with intercept and geno terms
      full<-control_em(initial = start[c(1:2,length(start))], phe = residuals2, desi = design[,1:2], pri = t(GP),
                       qua = qu, tole = tol, maxI = maxIter, count=count, logEM=logEM, W=W)

      ## calculating SE of effect size estimate for geno - for being able to get P value via wald test
      score_vectors <- calc_score_vector(esti = full$par, phen = residuals2, dsgn = design[,1:2], prior = t(GP), quan = qu, count = count)
      std_errors <- sqrt(diag(solve(t(score_vectors) %*% score_vectors)))

      full$chisq<-full$par[1]**2/(std_errors[1])**2
      full$Pval<-1-pchisq(full$chisq,df=1)
      full$df<-1

    } else{
      null<-control_em(initial = startNull, phe = pheno, desi = designNull, pri = t(GP),
                       qua = qu, tole = tol, maxI = maxIter, count=count, logEM=logEM, W=W)
      ## priming coeffecients of EM algorithm from null model
      start2<-c(start[1],null$par)
      
      full<-control_em(initial = start2, phe = pheno, desi = design, pri = t(GP),
                       qua = qu, tole = tol, maxI = maxIter, count=count, logEM=logEM, W=W)
    }

    if(SNPTEST){
      dfs<-full$df
      m2lq<-full$chisq
      pval<-full$Pval
      test <- c(df = dfs, m2lq = m2lq, pval = pval)        
    } else{
      dfs <- length(full$par) - length(null$par)
      m2lq <- 2*(null$value-full$value)
      pval <- pchisq(q = m2lq, df = dfs,lower.tail = FALSE)
      test <- c(df = dfs, m2lq = m2lq, pval = pval)
    }
    
    return_list <- list(full = full, test = test)
    if(se){ #assumes that first model is the full model
        score_vectors <- calc_score_vector(esti = full$par, phen = pheno, dsgn = design, prior = t(GP), quan = qu, count = count)
        std_errors <- sqrt(diag(solve(t(score_vectors) %*% score_vectors)))
        return_list$se <- std_errors
    }
    return(return_list)
}

#' Calculate minus log of likelihood used for dosage tests for quantitative traits
#'
#' @param params Parameters (first genetic effects, then intercept and additional covariate effects and last sd)
#' @param phenos Individual phenotypes length N
#' @param design Expected design matrix, N rows
#' @return Minus log of likelihood
dosage_mll_lm <- function(params, phenos, design) {
  probs<-dnorm(x = phenos,
               mean = design%*%(params[-length(params)]),
               sd = params[length(params)])
  return(-sum(log(probs)))
}

#' Calculate minus log of likelihood used for dosage tests for dichotomuous traits
#'
#' @param params Parameters (first genetic effects, then intercept and additional covariate effects)
#' @param phenos Individual phenotypes length N
#' @param design Expected design matrix N rows
#' @return Minus log of likelihood
dosage_mll_bin <- function(params, phenos, design) {
  probs <- dbinom(x=phenos,
                  size=1,
                  prob=plogis(design%*%params))
  return(-sum(log(probs)))
}

#' Calculate minus log of likelihood used for dosage tests for count traits
#'
#' @param params Parameters (first genetic effects, then intercept and additional covariate effects)
#' @param phenos Individual phenotypes length N
#' @param design Expected design matrix N rows
#' @return Minus log of likelihood
dosage_mll_pois <- function(params, phenos, design) {
  probs <- dpois(x=phenos,
                  lambda=exp(design%*%params))
  return(-sum(log(probs)))
}


#' Estimate parameters and calculate minus log of likelihood used for dosage tests
#'
#' @param params Parameters (first genetic effects, then intercept and additional covariate effects (and last sd for quantitative traits)
#' @param phenos Individual phenotypes length N
#' @param design Expected design matrix N rows
#' @return List of Parameters (par) and Minus log of likelihood (value)
dosage_fit <- function(phenos, design, qu, count){

    if(count){        
        fit <- fit_pois(design = design, pheno = phenos)
        para <- fit$coefficients
        mll <- dosage_mll_pois(params = para, phenos = phenos, design = design)        
    } else if(qu){     
        fit <- fit_lm(design = design, pheno = phenos)        
        sigm <- sqrt(sum(fit$residuals^2) / fit$df.residual)
        para <- c(as.vector(fit$coefficients), sigm)
        mll <- dosage_mll_lm(params = para, phenos = phenos, design = design)
    }else{
        fit <- fit_bin(design = design, pheno = phenos)
        para <- fit$coefficients
        mll <- dosage_mll_bin(params = para, phenos = phenos, design = design)
    }
    return(list(par = para, value = mll))
}

ANGSDasso_dosage <-function(GP, ys, # mandatory arguments 
                         covs, start, # can be missing
                         model='add', qu=TRUE, count=FALSE) {
  if(missing(covs)){
    covs<-matrix(1,nrow=length(ys),ncol=1)
  }
  
  if(ncol(GP)!=3 | nrow(GP)!=length(ys)){
    print("not genotype probs - has to have 3 cols and N rows")
    stop()
  }

  if(qu & count){
    print("have to specify either count or quantitative data")
    stop()
  }

  ## ASSUMES THAT INTERCEPT IS PART OF COVARIATES
  if(!all(covs[,1]==1)){
      print("seems like there might not be an intercept in covariates - HAS TO BE INCLUDED!!")
  }
   
  n<-length(ys)
  
  dosage<-colSums(t(GP)*c(0,1,2))
  
  if(model=='add'){    
    dosage<-colSums(t(GP)*c(0,1,2))
    design<-cbind(dosage,covs)
  } else if(model=='rec'){
    dosage<-colSums(t(GP)*c(0,0,1))
    design<-cbind(dosage,covs)
  } else{
    print("not valid model")
    stop()
  }
  
  if(missing(start)){
    start <- runif(ncol(design),-1,1)
    if(qu){
      start <- c(start,sd(ys))
    }
  }
  
  startNull<-start[-1]
  designNull<-as.matrix(design[,-1])
  
  null <- dosage_fit(phenos = ys, design = designNull, qu = qu, count=count)
  
  full <- dosage_fit(phenos = ys, design = design, qu = qu, count=count)
       
  dfs <- length(full$par) - length(null$par)
  m2lq <- 2*(null$value-full$value)
  pval <- pchisq(q = m2lq, df = dfs,lower.tail = FALSE)
  test <- c(df = dfs, m2lq = m2lq, pval = pval)
  
  return_list <- list(full = full, test = test)
  
  return(return_list)
}

