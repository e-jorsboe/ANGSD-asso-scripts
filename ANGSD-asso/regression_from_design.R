#' Fit linear regression model from design matrix in just R code
#'
#' @param design full design matrix
#' @param pheno responce variable
#' @param weights optional weights for weighted regression
#' @return list of coefficients, df.residuals and residuals (last two not meaningful for weights!?)
#' @examples
#' ## Annette Dobson (1990) "An Introduction to Generalized Linear Models".
#' ## Page 9: Plant Weight Data.
#' ctl <- c(4.17,5.58,5.18,6.11,4.50,4.61,5.17,4.53,5.33,5.14)
#' trt <- c(4.81,4.17,4.41,3.59,5.87,3.83,6.03,4.89,4.32,4.69)
#' group <- gl(2, 10, 20, labels = c("Ctl","Trt"))
#' weight <- c(ctl, trt)
#' lm_old <- lm(weight ~ group)
#' lm_new<-fit_lm(design = cbind(1, group=="Trt"), pheno = weight)
#' lm_new$coef
#' lm_old$coef
fit_lm <- function(design, pheno, weights=NULL){    
  if (dim(design)[1] != length(pheno)) { stop("dimention mismatch in fit_lm") }    
  if (is.null(weights)){
    weights <- rep(1, length(pheno))
  }
  z <- (weights == 0)
  xx <- design[!z, , drop=FALSE]
  yy <- pheno[!z]
  ww <- weights[!z]
  ttol <- 1e-07
  wwts <- sqrt(ww)
  #z <- .Call(stats:::C_Cdqrls, xx * wwts, yy * wwts, ttol, FALSE)
  xw <- xx * wwts
  yw<- yy * wwts  
  coef <- (solve(t(xw) %*% xw)) %*% (t(xw) %*% yw)   
  df <- length(yy) - length(coef)
  residuals <- pheno - design %*% coef
  return(list(coefficients = coef, df.residual = df, residuals = residuals))
}

#' Fit binomial regression model from design matrix in just R code
#'
#' @param design full design matrix
#' @param pheno responce variable
#' @param weights optional weights for weighted regression
#' @return list of coefficients, df.residuals and residuals (last two not meaningful for weights!?)

fit_bin <- function(design, pheno, weights=NULL){
  if (is.null(weights)) {
    weights <- rep(1, length(pheno))
  }
  z <- weights==0
  xx <- design[!z, , drop=FALSE]
  yy <- pheno[!z]
  ww <- weights[!z]
  ttol <- 1e-07
  
  linkinv <- function(x) {
    as.vector(1/(1 + exp(-(x))))
  }
  link <- function(x) {
    log(x/(1 - x))
  }
  variance <- function(mu){
    mu * (1 - mu)
  }
  coef <- rep(0, ncol(xx)) #estimates
  
  mustart <- (ww * yy + 0.5)/(ww + 1)
  eta <- link(mustart)
  mu <- linkinv(eta)
  
  mu.eta <- function(x){
    variance(linkinv(x))
  }
  for (i in 1:20) {
    varmu <- variance(mu)
    mu.eta.val <- mu.eta(eta)    
    z <- (eta) + (yy - mu)/mu.eta.val
    w <- sqrt((ww * mu.eta.val^2)/variance(mu))        
    ## fit <- .Call(stats:::C_Cdqrls, xx * w, z * w, 1e-07, check = FALSE)
    ## coef <- fit$coefficients
    xw <- xx * w
    yw <- z * w
    coef <- (solve(t(xw) %*% xw)) %*% (t(xw) %*% yw)
    ## this is the log OR
    eta <- as.vector(xx %*% coef)
    ## this is p(x)
    mu <- linkinv(eta)
  }

  #glm(pheno~.-1,data=dat, w=weights, family=binomial())$coefficients
  df <- length(yy) - length(coef)
  eta <- design %*% coef
  mu <- linkinv(eta)
  phenoW <- pheno
  phenoW[weights==0] <- 0
  residuals <- (phenoW - mu)/mu.eta(eta)
  return(list(coefficients = coef, df.residual = df, residuals = residuals))
}



#' Fit poisson regression model from design matrix in just R code
#'
#' @param design full design matrix
#' @param pheno responce variable
#' @param weights optional weights for weighted regression
#' @return list of coefficients, df.residuals and residuals (last two not meaningful for weights!?)

fit_pois <- function(design, pheno, weights=NULL){
  if (is.null(weights)) {
    weights <- rep(1, length(pheno))
  }
  z <- weights==0
  xx <- design[!z, , drop=FALSE]
  yy <- pheno[!z]
  ww <- weights[!z]
  ttol <- 1e-07
  ## changed this
  ## going from XB = g(u) to u = g^-1(XB)
  linkinv <- function(x) {
    as.vector(exp(x))
  }
  ## changed this, this is the link function
  ## or the 'g' function
  link <- function(x) {
      log(x)
  }
  ## changed this
  ## the variance in possion is given by lambda
  ## which is the same as our left side or mu = lambda
  variance <- function(mu){
      mu
  }
  coef <- rep(0, ncol(xx)) #estimates  
  mustart <- (ww * yy + 0.5)/(ww + 1)
  eta <- link(mustart)
  mu <- linkinv(eta)  
  mu.eta <- function(x){
      variance(linkinv(x))
  }
  for (i in 1:20) {
    varmu <- variance(mu)
    mu.eta.val <- mu.eta(eta)    
    z <- (eta) + (yy - mu)/mu.eta.val
    w <- sqrt((ww * mu.eta.val^2)/variance(mu))        
    ## fit <- .Call(stats:::C_Cdqrls, xx * w, z * w, 1e-07, check = FALSE)
    ## coef <- fit$coefficients
    xw <- xx * w
    yw <- z * w    
    coef <- (solve(t(xw) %*% xw)) %*% (t(xw) %*% yw)
    eta <- as.vector(xx %*% coef)
    mu <- linkinv(eta)
  }

  #glm(pheno~.-1,data=dat, w=weights, family=binomial())$coefficients
  df <- length(yy) - length(coef)
  eta <- design %*% coef
  mu <- linkinv(eta)
  phenoW <- pheno
  phenoW[weights==0] <- 0
  residuals <- (phenoW - mu)/mu.eta(eta)
  return(list(coefficients = coef, df.residual = df, residuals = residuals))
}
