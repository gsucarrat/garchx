####################################################
##
## CONTENTS:
##               
## TO DO LIST:
##
## ADD nobs and nonumber arguments to toLatex.garchx?
## ADD CHRISTOPHERSEN TEST?
## ADD more densities?
## ADD power unequal to 2?
##        
## 1 load libraries, compile C-code, create generics
##   GARCHRECURSION      #c-code recursion
##   refit               #generic for S3 methods
##
## 2 main garchx functions:
##   garchxSim           #simulate
##   garchxAvar          #compute asymptotic coefficient-covariance
##   garchxRecursion     #recursion
##   garchxObjective     #average log-likelihood
##   garchx              #estimate
##
## 3 extraction functions:
##   coef.garchx
##   fitted.garchx
##   logLik.garchx
##   nobs.garchx
##   predict.garchx
##   print.garchx
##   quantile.garchx
##   residuals.garchx
##   refit.garchx
##   toLatex.garchx
##   vcov.garchx
##
## 4 additional functions:
##   glag      #lag variable 
##   gdiff     #diff variable
##   refit     #generic for S3 methods
##   rmnorm    #draw from multivariate normal
##   ttest0    #t-test under nullity
##   waldtest0 #wald-test under nullity
##
####################################################


####################################################
## 1 LOAD LIBRARIES, STARTUP MESSAGE, COMPILE C-CODE
####################################################

##load required libraries
##=======================

library(zoo)

###C-code (for development-mode):
###===============================
#
#library(inline)
#
###garch recursion:
#mysig <- signature(iStart="integer", iEnd="integer",
#  iGARCHorder="integer", sigma2="numeric", parsgarch="numeric",
#  innov="numeric")
#mycode <- "
#
#  double garchsum; 
#  for(int i=*iStart; i < *iEnd; i++){
#
#    /*GARCH sum:*/
#    garchsum = 0;
#    for(int j=0; j < *iGARCHorder; j++){
#      garchsum = garchsum + parsgarch[j] * sigma2[i-1-j];
#    }
#
#    /*recursion:*/
#    sigma2[i] =  garchsum + innov[i];
#
#  }
#"
#GARCHXRECURSION <- cfunction(mysig, mycode, convention=".C") #compile
#rm("mysig", "mycode")
#
##garch recursion for simulations:
#mysig <- signature(iStart="integer", iEnd="integer",
#  iARCHorder="integer", iGARCHorder="integer", iASYMorder="integer",
#  parsarch="numeric", parsgarch="numeric", parsasym="numeric",
#  sigma2="numeric", z2="numeric", Ineg="numeric", xregsum="numeric")
#mycode <- "
#
#  double archsum;
#  double garchsum;
#  double asymsum;
#  
#  archsum = 0;
#  garchsum = 0;
#  asymsum = 0;
#
#  for(int i=*iStart; i < *iEnd; i++){
#
#    /* ARCH sum */
#    if(iARCHorder > 0){
#      archsum = 0;
#      for(int j=0; j < *iARCHorder; j++){
#        archsum = archsum + parsarch[j] * z2[i-1-j] * sigma2[i-1-j];
#      }
#    }
#
#    /* GARCH sum */
#    if(iGARCHorder > 0){
#      garchsum = 0;
#      for(int j=0; j < *iGARCHorder; j++){
#        garchsum = garchsum + parsgarch[j] * sigma2[i-1-j];
#      }
#    }
#    
#    /* ASYM sum */
#    if(iASYMorder > 0){
#      asymsum = 0;
#      for(int j=0; j < *iASYMorder; j++){
#        asymsum = asymsum + parsasym[j] * Ineg[i-1-j] * z2[i-1-j] * sigma2[i-1-j];
#      }
#    }
#    
#    /*recursion:*/
#    sigma2[i] = archsum + garchsum + asymsum + xregsum[i];
#	
#	} /* close for loop */
#
#}
#"
#GARCHXRECURSIONSIM <- cfunction(mysig, mycode, convention=".C") #compile
#rm("mysig", "mycode")

##create generics:
##================

##S3 generic/method 'refit':
refit <- function(object, ...){ UseMethod("refit") }


####################################################
## 2 MAIN FUNCTIONS
####################################################

##==================================================
## simulate from garch-x model:
garchxSim <- function(n, intercept=0.2, arch=0.1, garch=0.8, asym=NULL,
  xreg=NULL, innovations=NULL, backcast.values=list(), verbose=FALSE,
  as.zoo=TRUE, c.code=TRUE)
{
  ##orders:
  archOrder <- length(arch)
  garchOrder <- length(garch)
  asymOrder <- length(asym)
  maxOrder <- max(archOrder, garchOrder, asymOrder)

  ##initiate:
  if(is.null(innovations)){ innovations <- rnorm(n) }
  z2 <- innovations^2
  Ineg <- as.numeric(innovations < 0)
  sigma2 <- rep(0,n)
  if(is.null(xreg)){ xreg <- rep(0,n) }

  ##backcast values
  ##===============
  if(maxOrder > 0){

    ##innovations:
    if(is.null(backcast.values$innovations)){
      backcast.values$innovations <- rep(0, maxOrder)
    }
    innovations <- c(backcast.values$innovations, innovations)

    ##z2:
    if(is.null(backcast.values$z2)){
      z2mean <- mean(z2)
      backcast.values$z2 <- rep(z2mean, maxOrder)
    }
    z2 <- c(backcast.values$z2, z2)

    ##Ineg:
    if(is.null(backcast.values$Ineg)){
      backcast.values$Ineg <- rep(0, maxOrder)
    }
    Ineg <- c(backcast.values$Ineg, Ineg)

    ##sigma2:
    if(is.null(backcast.values$sigma2)){
      Esigma2 <- intercept/(1-sum(arch)-sum(garch))
      if( abs(Esigma2)==Inf ){ stop("Initial values of sigma2 are not finite") }
      backcast.values$sigma2 <- rep(Esigma2, maxOrder)
    }
    sigma2 <- c(backcast.values$sigma2, sigma2)

    ##xreg:
    if(is.null(backcast.values$xreg)){
      xregmean <- mean(xreg)
      backcast.values$xreg <- rep(xregmean, maxOrder)
    }
    xreg <- c(backcast.values$xreg, xreg)

  }

  ##recursion
  ##=========
  xregsum <- intercept + xreg

  if( c.code ){

    if( archOrder==0 ){ arch <- 0 }
    if( garchOrder==0 ){ garch <- 0 }
    if( asymOrder==0 ){ asym <- 0 }
    
#    ##development version (w/inline package):
#    tmp <- GARCHXRECURSIONSIM(as.integer(maxOrder), as.integer(length(sigma2)),
#      as.integer(archOrder), as.integer(garchOrder), as.integer(asymOrder),
#      as.numeric(arch), as.numeric(garch), as.numeric(asym),
#      as.numeric(sigma2), as.numeric(z2), as.numeric(Ineg),
#      as.numeric(xregsum))

    ##package version:
    tmp <- .C("GARCHXRECURSIONSIM", iStart=as.integer(maxOrder),
      iEnd=as.integer(length(sigma2)), iARCHorder=as.integer(archOrder),
      iGARCHorder=as.integer(garchOrder), iASYMorder=as.integer(asymOrder),
      parsarch=as.double(arch), parsgarch=as.double(garch),
      parsasym=as.double(asym), sigma2=as.double(sigma2),
      z2=as.double(z2), Ineg=as.double(Ineg), xregsum=as.double(xregsum),
      PACKAGE = "garchx")

    ##sigma2:
    sigma2 <- tmp$sigma2

  }else{

    archsum <- garchsum <- asymsum <- 0
    for(i in c(1+maxOrder):length(sigma2) ){
      if(archOrder > 0){
        archsum <-
          sum(arch*z2[c(i-1):c(i-archOrder)]*sigma2[c(i-1):c(i-archOrder)])
      }
      if(garchOrder > 0){
        garchsum <- sum(garch*sigma2[c(i-1):c(i-garchOrder)])
      }
      if(asymOrder > 0){
        asymsum <- sum(asym*Ineg[c(i-1):c(i-asymOrder)]*z2[c(i-1):c(i-asymOrder)]*sigma2[c(i-1):c(i-asymOrder)])
      }
      sigma2[i] <- archsum + garchsum + asymsum + xregsum[i]
    } #end for-loop

  } #end if( c.code )

  ##prepare result
  ##==============

  if(verbose){
    sigma <- sqrt(sigma2)
    y <- sigma*innovations
    result <- cbind(y, sigma, sigma2, Ineg, innovations)
    if(maxOrder > 0){ result <- result[-c(1:maxOrder),] }
    ##ensure result is a matrix when n=1:
    if(n==1){
      result <- rbind(result)
      rownames(result) <- NULL
    }
  }else{
    sigma <- sqrt(sigma2)
    result <- sigma*innovations
    if(maxOrder > 0){ result <- result[-c(1:maxOrder)] }
  }

  ##return result
  ##=============

  if(as.zoo){ result <- as.zoo(result) }
  return(result)

} ## close garchxSim function


##==================================================
## asymptotic coefficient-covariance
garchxAvar <- function(pars, arch=NULL, garch=NULL, asym=NULL,
  xreg=NULL, vcov.type=c("ordinary", "robust", "hac"),  innovations=NULL,
  Eeta4=NULL, n=1000000, objective.fun=1, seed=NULL)
{
  ##arguments
  ##---------

  ##pars:
  if( length(pars)==0 ){ stop("length(pars) cannot be 0") }

  ##vcov.type:
  types <- c("ordinary", "robust", "hac")
  whichType <- charmatch(vcov.type[1], types)
  vcov.type <- types[whichType]
  if( vcov.type %in% c("robust","hac") ){ stop("Sorry, not implemented yet!") }

  ##initiate
  ##--------
  parnames <- "intercept"
  aux <- list()
  aux$y.n <- n
  aux$recursion.n <- aux$y.n #used in recursion only; needed for robust vcov
  aux$y.coredata <- 1 #preliminary placeholder
  aux$y2 <- aux$y.coredata^2
  aux$y2mean <- 1 #preliminary placeholder

  ##arch, garch, asym arguments
  ##---------------------------

  ## note: the K refers to how many arch/asym/garch terms there are,
  ## NOT the arch/asym/garch order

  ##arch:
  if(is.null(arch)){
    aux$archK <- 0
  }else{
    if( identical(arch,0) ){ arch <- NULL }
    aux$archK <- length(arch)
  }
  aux$arch <- arch
  aux$archOrder <- ifelse(is.null(aux$arch), 0, max(aux$arch))

  ##garch:
  if(is.null(garch)){
    aux$garchK <- 0
  }else{
    if( identical(garch,0) ){ garch <- NULL }
    aux$garchK <- length(garch)
  }
  aux$garch <- garch
  aux$garchOrder <- ifelse(is.null(aux$garch), 0, max(aux$garch))

  ##asym:
  if(is.null(asym)){
    aux$asymK <- 0
  }else{
    if( identical(asym,0) ){ asym <- NULL }
    aux$asymK <- length(asym)
  }
  aux$asym <- asym
  aux$asymOrder <- ifelse(is.null(aux$asym), 0, max(aux$asym))

  ##xregK, maxpqr, maxpqrpluss1:
  aux$xregK <- ifelse(is.null(xreg), 0, NCOL(xreg))
  aux$maxpqr <- max(aux$archOrder,aux$garchOrder,aux$asymOrder)
  aux$maxpqrpluss1 <- aux$maxpqr + 1

  ##parameter indices and names
  ##---------------------------

  ##arch:
  if( aux$archK>0 ){
    aux$archIndx <- 2:c( length(arch)+1 )
    parnames <- c(parnames, paste0("arch", aux$arch))
  }else{ aux$archIndx <- 0 }

  ##garch:
  if( aux$garchK>0 ){
    aux$garchIndx <-
      c( max(1,aux$archIndx) +1) :
      c( max(1,aux$archIndx)+ length(aux$garch) )
    parnames <- c(parnames, paste0("garch", aux$garch))
  }else{ aux$garchIndx <- 0 }

  ##asym:
  if( aux$asymK>0 ){
    aux$asymIndx <-
      c( max(1,aux$archIndx,aux$garchIndx) +1) :
      c( max(1,aux$archIndx,aux$garchIndx) + length(aux$asym) )
    parnames <- c(parnames, paste0("asym", aux$asym))
  }else{ aux$asymIndx <- 0 }

  ##xreg:
  if( aux$xregK>0 ){
    aux$xregIndx <-
      c( max(1,aux$archIndx,aux$garchIndx,aux$asymIndx) +1) :
      c( max(1,aux$archIndx,aux$garchIndx,aux$asymIndx) +aux$xregK)
    xregNames <- colnames(xreg)
    if(is.null(xregNames)){
      xregNames <- paste0("xreg", 1:aux$xregK)
    }
    alreadyTaken <- c(1:aux$xregK)[ xregNames %in% parnames ]
    if( length(alreadyTaken)>0 ){
      xregNames[ alreadyTaken ] <- ""
    }
    missingNames <- which( xregNames=="" )
    if( length(missingNames) > 0 ){
      for(i in missingNames){
        xregNames[i] <- paste0("xreg", i)
      }
    }
    parnames <- c(parnames, xregNames)
  }else{ aux$xregIndx <- 0 }

  ##simulate
  ##--------

  ##arch argument:
  if( aux$archK==0 ){
    archArg <- NULL
  }else{
    archArg <- rep(0, aux$archOrder)
    archArg[ aux$arch ] <- pars[ aux$archIndx ]
  }

  ##garch argument:
  if( aux$garchK==0 ){
    garchArg <- NULL
  }else{
    garchArg <- rep(0, aux$garchOrder)
    garchArg[ aux$garch ] <- pars[ aux$garchIndx ]
  }

  ##asym argument:
  if( aux$archK==0 ){
    asymArg <- NULL
  }else{
    asymArg <- rep(0, aux$asymOrder)
    asymArg[ aux$asym ] <- pars[ aux$asymIndx ]
  }

  ##xreg argument:
  if( aux$xregK==0 ){
    xregArg <- NULL
  }else{
    xregArg <- cbind(xreg) %*% pars[ aux$xregIndx ]
  }

  ##set simulation length:
  if( aux$garchOrder>0 && is.null(innovations) ){
    nadj <- n+aux$garchOrder
  }else{ nadj <- n }

  ##simulate:
  if( !is.null(seed) ){ set.seed(seed) }
  mY <- garchxSim(nadj, intercept=pars[1], arch=archArg, garch=garchArg,
    asym=asymArg, xreg=xregArg, innovations=innovations, verbose=TRUE)
  if( nadj > n ){
    backcast.values <- as.numeric(mY[1:aux$garchOrder,"sigma2"])
    mY <- mY[ -c(1:aux$garchOrder),]
  }else{
    backcast.values <- NULL
  }

  ##y, y2, y2mean:
  aux$y.coredata <- coredata(mY[,"y"])
  aux$y2 <- aux$y.coredata^2
  aux$y2mean <- mean(aux$y2) #default garch backcast value

  ##auxiliary vectors and matrices
  ##------------------------------

  ##short y2, short ynotzero
  aux$y2short <- aux$y2[ aux$maxpqrpluss1:aux$y.n ]
  if( objective.fun==0 ){
    aux$ynotzero <-
      as.numeric(aux$y.coredata != 0)[ aux$maxpqrpluss1:aux$y.n ]
  }

  ##arch matrix:
  if( aux$archK>0 ){
    aux$y2matrix <- matrix(aux$y2mean, aux$y.n, aux$archK)
    for(i in 1:aux$archK){
      aux$y2matrix[ c(1+aux$arch[i]):aux$y.n ,i] <-
        aux$y2[ 1:c(aux$y.n-aux$arch[i]) ]
    }
  }

  ##garch vector:
  if( aux$garchK>0 ){
    aux$sigma2 <- rep(aux$y2mean, aux$y.n)
    if( !is.null(backcast.values) ){
      aux$sigma2[1:aux$garchOrder] <- backcast.values
    }
  }

  ##asym matrix:
  if( aux$asymK>0 ){
    aux$Ineg <- as.numeric(aux$y.coredata < 0)
    aux$Inegy2 <- aux$Ineg*aux$y2
    backvals <- mean(aux$Inegy2)
    aux$Inegy2matrix <- matrix(backvals, aux$y.n, aux$asymK)
    for(i in 1:aux$asymK){
      aux$Inegy2matrix[ c(1+aux$asym[i]):aux$y.n ,i] <-
        aux$Inegy2[ 1:c(aux$y.n-aux$asym[i]) ]
    }
  }

  ##xreg:
  if(aux$xregK>0){
    colnames(xreg) <- NULL
    aux$xreg <- cbind(xreg)
  }

  ##miscellaneous
  ##-------------

  aux$upper <- Inf
  aux$lower <- 0
  aux$control <- list()
  aux$hessian.control <- list()
  aux$solve.tol <- .Machine$double.eps
  aux$c.code <- TRUE
  aux$sigma2.min <- .Machine$double.eps
  aux$objective.fun <- objective.fun
  aux$penalty.value <- garchxObjective(pars, aux)

  ##hessian
  ##-------

  names(pars) <- parnames
  aux$hessian <- optimHess(pars, garchxObjective, aux=aux,
    control=aux$hessian.control)

  ##ordinary vcov
  ##-------------
  if( vcov.type=="ordinary" ){

    ##kappa:
    if( is.null(Eeta4) ){
      if( aux$objective.fun==1 ){
        Eeta4 <- mean( mY[,"innovations"]^4 )
      }
      if( aux$objective.fun==0 ){
        etanotzero <- as.numeric( mY[,"innovations"] != 0 )
        Eeta4 <- sum( etanotzero*innovations^4 )/sum(etanotzero)
      }
    }

    ##avar:
    result <- (Eeta4 - 1)*solve(aux$hessian, tol=aux$solve.tol)
  }


  ##robust vcov
  ##-----------
  if( vcov.type=="robust" ){
    ##tba
  }

  ##hac vcov
  ##-----------
  if( vcov.type=="hac" ){
    ##tba
  }

  ##return
  ##------
  return(result)

} #close


##==================================================
## recursion of garch-model
garchxRecursion <- function(pars, aux)
{
  ##initiate/intercept:
  innov <- rep(pars[1], aux$y.n)
  
  ##arch:
  if( aux$archK>0 ){
    ##use crossprod() instead of %*%?
    innov <- innov + aux$y2matrix %*% pars[ aux$archIndx ]
  }

  ##asym:
  if( aux$asymK>0 ){
    ##use crossprod() instead of %*%?
    innov <- innov + aux$Inegy2matrix %*% pars[ aux$asymIndx ]
  }

  ##xreg:
  if( aux$xregK>0 ){
    ##use crossprod() instead of %*%?
    innov <- innov + aux$xreg %*% pars[ aux$xregIndx ]
  }
    

  ##garchK=0
  ##--------
  if( aux$garchK == 0 ){
    sigma2 <- as.vector(innov)
  }


  ##garchK>0
  ##--------  
  if( aux$garchK > 0 ){
  
    innov <- as.vector(innov)
    sigma2 <- aux$sigma2
    parsgarch <- rep(0, aux$garchOrder)
    parsgarch[ aux$garch ] <- pars[ aux$garchIndx ]
        
    ##if(c.code):
    if(aux$c.code){

#      ##development version (w/inline package):
#      tmp <- GARCHXRECURSION(as.integer(aux$garchOrder),
#        as.integer(aux$recursion.n), as.integer(aux$garchOrder),
#        as.numeric(sigma2), as.numeric(parsgarch),
#        as.numeric(innov))

      ##package version:
      tmp <- .C("GARCHXRECURSION", iStart = as.integer(aux$garchOrder),
        iEnd = as.integer(aux$recursion.n), iGARCHorder = as.integer(aux$garchOrder),
        sigma2 = as.double(sigma2), parsgarch = as.double(parsgarch),
        innov = as.double(innov), PACKAGE = "garchx")

      ##sigma2:
      sigma2 <- tmp$sigma2

    }else{

      ##if(garch1):
      if(aux$garchOrder==1){
        for( i in 2:aux$y.n ){
          sigma2[i] <- parsgarch * sigma2[i-1] + innov[i]
        }
      }else{
        for( i in c(1+aux$garchOrder):aux$recursion.n ){
          sigma2[i] <-
            sum(parsgarch*sigma2[ c(i-1):c(i-aux$garchOrder) ]) + innov[i]
        }
      }

    } #close if(c.code)

  } #close if(garchK==0)else(..)
  

  ##output:
  return(sigma2)
  ##note: this volatility contains the first observations
  ##(i.e. 1:aux$maxpqr), which should be deleted before
  ##returned by fitted(..) etc.

} ##close garchxRecursion function
  

##==================================================
## average loglikelihood
garchxObjective <- function(pars, aux)
{
  ##initiate:
  parsOK <- TRUE; sigma2OK <- TRUE
  
  ##check parameters for NAs:
  if( any(is.na(pars))){ parsOK <- FALSE  }
  
  ##check xreg parameters:
  if( aux$xregK>0 ){
    value <- pars[1] + sum(pars[aux$xregIndx])
    if( value <= 0 ){ parsOK <- FALSE }
  }

  ##compute and check sigma2:
  if( parsOK ){
    sigma2 <- garchxRecursion(pars, aux)
    sigma2 <- sigma2[ aux$maxpqrpluss1:aux$y.n ]
    if( any(sigma2<=0) || any(sigma2==Inf) ){
      sigma2OK <- FALSE
    }else{
      ##robustify log-transform in log-likelihood:
      sigma2 <- pmax( aux$sigma2.min, sigma2 )
    }
  }

  ##compute average log-likelihood:
  if( parsOK && sigma2OK ){
    if(aux$objective.fun==0){
      result <- mean( aux$ynotzero*(aux$y2short/sigma2 + log(sigma2)) )
    }
    if(aux$objective.fun==1){
      result <- mean( aux$y2short/sigma2 + log(sigma2) )
    }  
  }else{
    result <- aux$penalty.value
  }

  ##return result:
  return(result)

} ##close garchxObjective function

    
##==================================================
## estimate garch
garchx <- function(y, order=c(1,1), arch=NULL, garch=NULL, asym=NULL,
  xreg=NULL, vcov.type=c("ordinary","robust","hac"), initial.values=NULL,
  backcast.values=NULL, lower=0, upper=+Inf, control=list(),
  hessian.control=list(), solve.tol=.Machine$double.eps, estimate=TRUE,
  c.code=TRUE, penalty.value=NULL, sigma2.min=.Machine$double.eps,
  objective.fun=1, turbo=FALSE)
{
  ##sys.call:
  sysCall <- sys.call()
    
  ##create auxiliary list, date, parnames:
  aux <- list()
  aux$date <- date()
  aux$sys.call <- sysCall
  parnames <- "intercept"

  ##y argument
  ##----------
  
  aux$y.name <- deparse(substitute(y))
  y <- na.trim(as.zoo(y))
  aux$y.n <- NROW(y)
  aux$recursion.n <- aux$y.n #used in recursion only; needed for robust vcov
  aux$y.coredata <- as.vector(coredata(y)) #in case y is matrix (e.g. due to xts)
  aux$y.index <- index(y)
  aux$y2 <- aux$y.coredata^2
  aux$y2mean <- mean(aux$y2) #default garch backcast value if is.null(backcast.values)

  ##order argument
  ##--------------

  aux$order <- c(0,0,0)
  if( length(order)>0 ){
    aux$order[ 1:length(order) ] <- order[ 1:length(order) ]
  }
  if(is.null(garch) && aux$order[1]>0){ garch <- 1:aux$order[1] }
  if(is.null(arch) && aux$order[2]>0){ arch <- 1:aux$order[2] }
  if(is.null(asym) && aux$order[3]>0){ asym <- 1:aux$order[3] }
  
  ##arch, garch, asym arguments
  ##---------------------------

  ## note: the K refers to how many arch/asym/garch terms there are,
  ## NOT the arch/asym/garch order

  ##arch:
  if(is.null(arch)){
    aux$archK <- 0
  }else{
    if( identical(arch,0) ){ arch <- NULL }  
    aux$archK <- length(arch)
  }
  aux$arch <- arch
  aux$archOrder <- ifelse(is.null(aux$arch), 0, max(aux$arch))

  ##garch:
  if(is.null(garch)){
    aux$garchK <- 0
  }else{
    if( identical(garch,0) ){ garch <- NULL }  
    aux$garchK <- length(garch)
  }
  aux$garch <- garch
  aux$garchOrder <- ifelse(is.null(aux$garch), 0, max(aux$garch))
    
  ##asym:
  if(is.null(asym)){
    aux$asymK <- 0
  }else{
    if( identical(asym,0) ){ asym <- NULL }  
    aux$asymK <- length(asym)
  }
  aux$asym <- asym
  aux$asymOrder <- ifelse(is.null(aux$asym), 0, max(aux$asym))

  ##xregK, maxpqr, maxpqrpluss1:
  aux$xregK <- ifelse(is.null(xreg), 0, NCOL(xreg))
  aux$maxpqr <- max(aux$archOrder,aux$garchOrder,aux$asymOrder)
  aux$maxpqrpluss1 <- aux$maxpqr + 1
  
  ##parameter indices and names
  ##---------------------------
  
  ##arch:
  if( aux$archK>0 ){
    aux$archIndx <- 2:c( length(arch)+1 )
    parnames <- c(parnames, paste0("arch", aux$arch))
  }else{ aux$archIndx <- 0 }
  
  ##garch:
  if( aux$garchK>0 ){
    aux$garchIndx <-
      c( max(1,aux$archIndx) +1) :
      c( max(1,aux$archIndx)+ length(aux$garch) )
    parnames <- c(parnames, paste0("garch", aux$garch))
  }else{ aux$garchIndx <- 0 }
  
  ##asym:
  if( aux$asymK>0 ){
    aux$asymIndx <-
      c( max(1,aux$archIndx,aux$garchIndx) +1) : 
      c( max(1,aux$archIndx,aux$garchIndx) + length(aux$asym) )
    parnames <- c(parnames, paste0("asym", aux$asym))
  }else{ aux$asymIndx <- 0 }
  
  ##xreg:
  if( aux$xregK>0 ){
    aux$xregIndx <-
      c( max(1,aux$archIndx,aux$garchIndx,aux$asymIndx) +1) : 
      c( max(1,aux$archIndx,aux$garchIndx,aux$asymIndx) +aux$xregK)
    xregNames <- colnames(xreg)
    if(is.null(xregNames)){
      xregNames <- paste0("xreg", 1:aux$xregK)
    }
    alreadyTaken <- c(1:aux$xregK)[ xregNames %in% parnames ]
    if( length(alreadyTaken)>0 ){
      xregNames[ alreadyTaken ] <- ""
    }
    missingNames <- which( xregNames=="" )
    if( length(missingNames) > 0 ){
      for(i in missingNames){
        xregNames[i] <- paste0("xreg", i)
      }
    }
    parnames <- c(parnames, xregNames)
  }else{ aux$xregIndx <- 0 }

  ##backcast value:
  ##---------------
  if( is.null(backcast.values) ){
    backcast.values <- aux$y2mean
  }    
  if( length(backcast.values)!=1 ) stop("length(backcast.values) must be 1")
  if( backcast.values < 0 ) stop("'backcast.values' must be non-negative")
  aux$backcast.values <- backcast.values
  
  ##auxiliary vectors and matrices
  ##------------------------------

  ##short y2, short ynotzero
  aux$y2short <- aux$y2[ aux$maxpqrpluss1:aux$y.n ]
  if(objective.fun==0){
    aux$ynotzero <-
      as.numeric(aux$y.coredata != 0)[ aux$maxpqrpluss1:aux$y.n ]
  }

  ##arch matrix:
  if( aux$archK>0 ){
    aux$y2matrix <- matrix(aux$y2mean, aux$y.n, aux$archK)
    for(i in 1:aux$archK){
      aux$y2matrix[ c(1+aux$arch[i]):aux$y.n ,i] <-
        aux$y2[ 1:c(aux$y.n-aux$arch[i]) ] 
    }
  }

  ##garch vector:
  if( aux$garchK>0 ){
    aux$sigma2 <- rep(aux$y2mean, aux$y.n)
    if( !is.null(backcast.values) ){
      aux$sigma2[1:aux$garchOrder] <- backcast.values
    }        
  }

  ##asym matrix:    
  if( aux$asymK>0 ){
    aux$Ineg <- as.numeric(aux$y.coredata < 0)
    aux$Inegy2 <- aux$Ineg*aux$y2
    backvals <- mean(aux$Inegy2)
    aux$Inegy2matrix <- matrix(backvals, aux$y.n, aux$asymK)
    for(i in 1:aux$asymK){
      aux$Inegy2matrix[ c(1+aux$asym[i]):aux$y.n ,i] <-
        aux$Inegy2[ 1:c(aux$y.n-aux$asym[i]) ] 
    }
  }

  ##xreg:  
  if(aux$xregK>0){
    xreg <- na.trim(as.zoo(xreg))
    xreg <-
      window(xreg, start=aux$y.index[1], end=aux$y.index[aux$y.n])
    xreg <- coredata(xreg)
    colnames(xreg) <- NULL
    aux$xreg <- cbind(xreg)
  }

 
  ##initial parameter values:
  ##-------------------------

  if(is.null(initial.values)){

    ##intercept:
    aux$initial.values <- 0.1

    ##arch:
    if( aux$archK>0 ){
      aux$initial.values <-
        c(aux$initial.values, rep(0.1/aux$archK, aux$archK)/aux$arch)
    }

    ##garch:
    if( aux$garchK>0 ){
      aux$initial.values <-
        c(aux$initial.values, rep(0.7/aux$garchK, aux$garchK)/aux$garch)
    }

    ##asym:
    if( aux$asymK>0 ){
      aux$initial.values <-
        c(aux$initial.values, rep(0.02, aux$asymK)/aux$asym)
    }

    ##xreg:
    if( aux$xregK>0 ){
      aux$initial.values <-
        c(aux$initial.values, rep(0.01, aux$xregK))
    }

  }else{
 
    aux$initial.values <- initial.values
    
  }
  
  ##here, I could add a check for stability: if unstable,
  ##stop or modify?


  ##miscellaneous
  ##-------------

  aux$upper <- upper
  aux$lower <- lower  
  aux$control <- control
  aux$hessian.control <- hessian.control
  aux$solve.tol <- solve.tol
  aux$c.code <- c.code
  aux$sigma2.min <- sigma2.min
  aux$objective.fun <- objective.fun
  aux$penalty.value <- penalty.value
  if(is.null(aux$penalty.value)){
    aux$penalty.value <- garchxObjective(aux$initial.values, aux)
  }

  ##estimation
  ##----------

  if(estimate){
    result <- nlminb(aux$initial.values, garchxObjective,
      aux=aux, control=aux$control, upper=aux$upper, lower=aux$lower)
  }else{
    result <- list()
    result$par <- aux$initial.values
    result$objective <- aux$penalty.value
    result$convergence <- NA
    result$iterations <- NA
    result$evaluations <- NA
    result$message <- "none, since 'estimate = FALSE'"
  }  
  names(result$par) <- parnames
  aux <- c(aux, result)
  names(aux$initial.values) <- parnames          

  ##not turbo?
  ##----------

  if( !turbo ){                               

    ##fitted values, residuals:
    sigma2 <- garchxRecursion(as.numeric(aux$par), aux)
    residStd <- aux$y.coredata/sqrt(sigma2)
    ##convert to zoo:
    sigma2 <- zoo(sigma2, order.by=aux$y.index)
    residStd <- zoo(residStd, order.by=aux$y.index)
    ##shorten vectors, add to result:
    aux$fitted <- sigma2[ aux$maxpqrpluss1:aux$y.n ]
    aux$residuals <- residStd[ aux$maxpqrpluss1:aux$y.n ]    

    ##hessian:    
    aux$hessian <- optimHess(aux$par, garchxObjective,
      aux=aux, control=aux$hessian.control)  

    ##vcov:
    aux$vcov <- vcov.garchx(aux, vcov.type=vcov.type)

  } #close if(!turbo)

  ##result
  ##------

  class(aux) <- "garchx"
  return(aux)

} ##close garchx() function


####################################################
## 3 EXTRACTION FUNCTIONS
####################################################

##==================================================
## extract coefficients
coef.garchx <- function(object, ...){ return(object$par) }

##==================================================
## compute confidence intervals of parameters
confint.garchx <- function(object, parm, level = 0.95, ...)
{
  ##the code of this function is based on and almost identical
  ##to that of confint.default()

  cf <- coef(object)
  pnames <- names(cf)
  if (missing(parm)) 
      parm <- pnames
  else if (is.numeric(parm)) 
      parm <- pnames[parm]
  a <- (1 - level)/2
  a <- c(a, 1 - a)
  
  ##lower and upper names:
  lowerName <- paste0(a[1]*100, " %")
  upperName <- paste0(a[2]*100, " %")
  pct <- c(lowerName, upperName)
  
  ##compute limits:
  iDF <- nobs.garchx(object) - length(coef.garchx(object))
  fac <- qt(a, df=iDF)
  ci <- array(NA, dim = c(length(parm), 2L), dimnames = list(parm, 
      pct))
  ses <- sqrt(diag(vcov(object)))[parm]
  ci[] <- cf[parm] + ses %o% fac
  
  ##return result:
  return(ci)

} #close confint.garchx()

##==================================================
## extract fitted values
fitted.garchx <- function(object, as.zoo=TRUE, ...){
  if(is.null(object$fitted)){
    object$fitted <- garchxRecursion(as.numeric(object$par), object)
    if(as.zoo){
      object$fitted <- zoo(object$fitted, order.by=object$y.index)
    }
    object$fitted <- object$fitted[ object$maxpqrpluss1:object$y.n ]
  }
  return(object$fitted)
}

##==================================================
## gaussian log-likelihood:
logLik.garchx <- function(object, ...){
  y <- object$y.coredata[ object$maxpqrpluss1:object$y.n ]
  sigma2 <- coredata(fitted.garchx(object))
  dnormvals <- dnorm(y, mean=0, sd=sqrt(sigma2), log=TRUE)
  if(object$objective.fun==0){
    result <- sum( object$ynotzero * dnormvals )
  }
  if(object$objective.fun==1){
    result <- sum(dnormvals)
  }
  attr(result, "df") <- length(object$par)
  attr(result, "nobs") <- ifelse(object$objective.fun==1,
    length(sigma2), sum(object$ynotzero))
  return(result)
}

##==================================================
## number of observations used in the estimation:
nobs.garchx <- function(object, ...){
  length(residuals.garchx(object))
}

##==================================================
## predict:
predict.garchx <- function(object, n.ahead=10, newxreg=NULL,
  newindex=NULL, n.sim=NULL, verbose=FALSE, ...)
{
  ##coefficients
  coefs <- as.numeric(coef.garchx(object))
  interceptCoef <- coefs[1]

  archCoef <- NULL
  if(object$archK>0){
    archCoef <- rep(0, object$archOrder)
    archCoef[ object$arch ] <- coefs[ object$archIndx ]
  }

  garchCoef <- NULL
  if(object$garchK>0){
    garchCoef <- rep(0, object$garchOrder)
    garchCoef[ object$garch ] <- coefs[ object$garchIndx ]
  }

  asymCoef <- NULL
  if(object$asymK>0){
    asymCoef <- rep(0, object$asymOrder)
    asymCoef[ object$asym ] <- coefs[ object$asymIndx ]
  }

  if( object$xregK>0 ){ xregCoef <- coefs[ object$xregIndx ] }


  ##backcast values
  backcast.values <- list()
  if( object$maxpqr>0 ){
    backcast.values$innovations <-
      tail(as.numeric(residuals.garchx(object)), n=object$maxpqr)
    backcast.values$z2 <- backcast.values$innovations^2
    backcast.values$Ineg <- as.numeric( backcast.values$innovations<0 )
    backcast.values$sigma2 <-
      tail(as.numeric(fitted.garchx(object)), n=object$maxpqr)
    if(object$xregK>0){
      backcast.values$xreg <-
        tail(as.numeric(object$xreg %*% xregCoef), n=object$maxpqr)
    }
  }

  ##newxreg
  xreg <- NULL
  if( object$xregK>0 ){
    if( (NROW(newxreg)!=n.ahead) ){
      stop("NROW(newxreg) is unequal to n.ahead")
    }
    xreg <- cbind(newxreg) %*% xregCoef
  }

  ##n.sim:
  if(is.null(n.sim)){
    if(n.ahead==1){ n.sim <- 1 }else{ n.sim <- 5000 }
  }

  ##bootstrap the innovations
  etahat <- coredata(residuals.garchx(object))
  draws <- runif(n.ahead * n.sim, min = 0.5 + 1e-09,
    max = length(etahat) + 0.5 - 1e-09)
  draws <- round(draws, digits = 0)
  innovations <- matrix(etahat[draws], n.ahead, n.sim)

  ##make predictions
  predictions <- matrix(NA, n.ahead, n.sim)
  for(i in 1:n.sim){
    predictions[,i] <- garchxSim(n.ahead, intercept=interceptCoef,
      arch=archCoef, garch=garchCoef, asym=asymCoef, xreg=xreg,
      innovations=innovations[,i], backcast.values=backcast.values,
      verbose=TRUE, as.zoo=FALSE)[,"sigma2"]
  }

  ##result
  result <- as.vector(rowMeans(predictions))
  if(verbose){
    result <- cbind(result, predictions)
    colnames(result) <- c("sigma2", 1:NCOL(predictions))
  }
  if(is.null(newindex)){ newindex <- 1:n.ahead }
  result <- zoo(result, order.by=newindex)
  return(result)

} #close predict.garchx


##==================================================
## print result:
print.garchx <- function(x, ...){
  
  ##out1:
  pars <- coef.garchx(x)
  vcovmat <- vcov.garchx(x)
  vcovComment <- comment(vcovmat)
  out1 <- rbind(pars, sqrt(diag(vcovmat)))
  rownames(out1) <- c("Estimate:", "Std. Error:")

  ##out2:
  out2 <- as.data.frame(matrix(NA, 1, 1))
  out2[1, 1] <- as.character(round(logLik.garchx(x), digits = 4))
  rownames(out2) <- "Log-likelihood:"
  colnames(out2) <- ""
    
  ##print message:
  cat("\n")
  cat("Date:", x$date, "\n")
#  cat("Dependent variable:", x$y.name, "\n")
  cat("Method: normal ML\n")
  cat("Coefficient covariance:", vcovComment, "\n")
  cat("Message (nlminb):", x$message, "\n")
  cat("No. of observations (fitted):", x$y.n - x$maxpqr, "\n")
  cat("Sample:", as.character(x$y.index[1]), "to", as.character(x$y.index[x$y.n]), 
      "\n")
  cat("\n")
#  cat("Coefficients: \n ")
  print(out1)
  print(out2)
  cat("\n")
}

##==================================================
## extract residuals:
quantile.garchx <- function(x, probs=0.025, names=TRUE,
  type=7, as.zoo=TRUE, ...)
{
  etahat <- residuals.garchx(x)
  sigma <- sqrt(fitted.garchx(x)) 
  qvals <- quantile(etahat, probs=probs, names=names, type=type) 
  iN <- NROW(etahat)
  iCols <- length(qvals)
  result <- matrix(NA, iN, iCols)
  colnames(result) <- names(qvals)
  for(i in 1:iCols){
    result[,i] <- sigma*qvals[i]
  }
  if(iCols==1){ result <- as.vector(result) }
  if(as.zoo){ result <- zoo(result, order.by=index(etahat)) }
  return(result)
}

##==================================================
## refit to new data:
refit.garchx <- function(object, newy=NULL, newxreg=NULL,
  backcast.value=NULL, reestimate=FALSE, ...)
{
  ##check:
  if( is.null(newy)){ stop("'newy' cannot be NULL") }
  if( !is.null(object$xreg) && is.null(newxreg) ){
    stop("'newxreg' is missing")
  }
  
  ##obtain garch spec:    
  archArg <- object$arch
  garchArg <- object$garch 
  asymArg <- object$asym

  ##refit:
  if(reestimate){
    result <- garchx(newy, arch=archArg, garch=garchArg, asym=asymArg,
      xreg=newxreg, backcast.values=backcast.value, ...)
  }else{
    coefs <- coef.garchx(object)
    result <- garchx(newy, arch=archArg, garch=garchArg, asym=asymArg,
      xreg=newxreg, initial.values=coefs, backcast.values=backcast.value,
      estimate=FALSE, turbo=TRUE)
    result$convergence <- NA
    result$iterations <- NA
    result$evaluations <- NA
    result$message <- "not applicable, since 'reestimate = FALSE'"
    result$fitted <- fitted.garchx(result)
    result$residuals <- residuals.garchx(result)
    result$hessian <- object$hessian
    result$vcov <- object$vcov
  }
    
  ##return result:
  return(result)
}

##==================================================
## extract residuals:
residuals.garchx <- function(object, as.zoo=TRUE, ...){
  if(is.null(object$residuals)){
    sigma2 <- garchxRecursion(as.numeric(object$par), object)
    object$residuals <- object$y.coredata/sqrt(sigma2)
    if(as.zoo){
      object$residuals <- zoo(object$residuals, order.by=object$y.index)
    }
    object$residuals <-
      object$residuals[ object$maxpqrpluss1:object$y.n ]
  }
  return(object$residuals)
}

##==================================================
## print LaTeX code (equation form):
toLatex.garchx <- function(object, digits=4, ...)
{
  ##equation:
  ##---------

  coefs <- coef.garchx(object)
  coefsNames <- names(coefs)
  coefsNames[1] <- "" #intercept
  coefs <- as.numeric(coefs)
  stderrs <- as.numeric(sqrt(diag(vcov(object))))

  eqtxt <- NULL
  for(i in 1:length(coefs) ){
    ifpluss <- ifelse(i==1, "", " + ")
    eqtxt <- paste(eqtxt,
      ifelse(coefs[i] < 0, " - ", ifpluss), "\\underset{(",
      format(round(stderrs[i], digits=digits), nsmall=digits), ")}{",
      format(round(abs(coefs[i]), digits=digits), nsmall=digits), "}",
      coefsNames[i], sep="")
  }

  txtAddEq <- " \\\\[1mm]"
  eqtxt <- paste0("  \\widehat{\\sigma}_t^2 &=& ", eqtxt, "",
    txtAddEq, " \n")

  ##goodness of fit:
  ##----------------

  goftxt <- NULL
  goftxt <- "   &&"
  iT <- nobs.garchx(object)
  goftxt <- paste(goftxt, " \\text{Log-likelihood: }",
    format(round(as.numeric(logLik.garchx(object)), digits=digits), nsmall=digits),
      "\\qquad T = ", iT, " \\nonumber \n", sep="")

  ##print code:
  ##-----------

  cat("%% the model was estimated", object$date, "\n")
  cat("%% note: the 'eqnarray' environment requires the 'amsmath' package\n")
  cat("\\begin{eqnarray}\n")
  cat(eqtxt)
  cat(goftxt)
  cat("\\end{eqnarray}\n")

} #close toLatex


##==================================================
## extract variance-covariance matrix:
vcov.garchx <- function(object, vcov.type=NULL, ...)
{
  ##compute hessian?
  ##----------------
  if( ! ("hessian" %in% names(object)) ){
    object$hessian <- optimHess(object$par, garchxObjective,
      aux=object, control=object$hessian.control)  
  }

  ##determine vcov.type
  ##-------------------
  vcovTypes <- c("ordinary", "robust", "hac")
  if( is.null(vcov.type) ){
    sysCall <- as.list(object$sys.call)
    if( "vcov.type" %in% names(sysCall) ){
      whichOne <- which( "vcov.type" == names(sysCall) )
      vcov.type <- sysCall[[ whichOne ]]
    }else{
      vcov.type <- vcovTypes
    }
  }    
  whichType <- charmatch(vcov.type[1], vcovTypes)
  vcov.type <- vcovTypes[ whichType ]

  ##vcov comment
  ##------------
  vcovComment <-
    ifelse("vcov" %in% names(object), comment(object$vcov), "NULL")

  ##ordinary vcov
  ##-------------
  if( vcov.type=="ordinary" && vcovComment != "ordinary" ){

    ##kappahat:
    if( is.null(object$residuals) ){
      object$residuals <- residuals.garchx(object)
    }
    if( object$objective.fun==0 ){
      kappahat <-
        sum( object$ynotzero * object$residuals^4)/sum(object$ynotzero)
    }
    if( object$objective.fun==1 ){
      kappahat <- mean( object$residuals^4 )
    }

    ##vcov:
    iN <- length(object$residuals) #divide by n for finite sample version
    object$vcov <-
      (kappahat - 1)*solve(object$hessian, tol=object$solve.tol)/iN
    comment(object$vcov) <- "ordinary"

  } #close if(ordinary)
        
  ##"robust" vcov
  ##-------------
  if( vcov.type=="robust" && vcovComment!="robust" ){

    ##inverse of J:
    Jinv <- solve(object$hessian, tol=object$solve.tol)
    
    ##temporary function:
    funtmp <- function(i, pars){
      object$recursion.n <- i
      object$par <- pars
      sigma2 <- garchxRecursion(as.numeric(object$par), object)
      return(sigma2[i])
    }

    ##matrix w/gradients
    pars <- as.numeric(object$par)
    mIadj <- matrix(NA, object$y.n, length(object$par))
#    sigma2Check <- rep(NA, object$y.n)
    for(i in 1:object$y.n){
      mIadj[i,] <-
        attr(numericDeriv(quote(funtmp(i,pars)),"pars"), "gradient")
#      ##code that enables a of sigma2:
#      tmp <- numericDeriv(quote(funtmp(i,pars)),"pars")
#      sigma2Check[i] <- tmp[1]
#      mIadj[i,] <- attr(tmp, "gradient")
    }

    ##Iadj:
    sigma2 <- garchxRecursion(as.numeric(object$par), object)
    etahatadj <- object$y2/(sigma2^2)
    mIadj <- etahatadj * mIadj
    mIadj <- mIadj[ object$maxpqrpluss1:object$y.n, ]
    iN <- length(residuals.garchx(object))    
    mIadj <- crossprod(mIadj)/iN - object$hessian

    ##vcov:
    object$vcov <- (Jinv %*% mIadj %*% Jinv)/iN
    comment(object$vcov) <- "robust"

  } #close if("robust")
    
  ##"hac" vcov
  ##----------
  if( vcov.type=="hac" && vcovComment!="hac" ){

    ##inverse of J:
    Jinv <- solve(object$hessian, tol=object$solve.tol)
    
    ##temporary function:
    ##(to be used to compute the derivatives of sigma2 at i)
    funtmp <- function(i, pars){
      object$recursion.n <- i
      object$par <- pars
      sigma2 <- garchxRecursion(as.numeric(object$par), object)
      return(sigma2[i])
    }

    ##make mldot: matrix w/gradients of sigma2
    pars <- as.numeric(object$par)
    mldot <- matrix(NA, object$y.n, length(object$par))
    for(i in 1:object$y.n){
      mldot[i,] <-
        attr(numericDeriv(quote(funtmp(i,pars)),"pars"), "gradient")
    }

    ##make mShat: matrix w/gradients of lhat
    sigma2 <- garchxRecursion(as.numeric(object$par), object)
    etahatadj <- object$y2/(sigma2^2)
    mShat <- (1/sigma2 - etahatadj) * mldot #matrix w/gradients of lhat

    ##bandwidth:
    igamma_n <- 4*(object$y.n/100)^(2/9) #bandwidth
    ##maximum lag:
    iL <- floor(igamma_n) #note: j/iS_T > 1 when j > i
    ##Bartlett weights:
    vW <- 1 - 1:iL/igamma_n #Hansen (1992)

    ##compute mJhat:
    iT <- NROW(mShat) #for more compact notation
    mIhat <- crossprod(mShat)/iT  #(j=0)
    if( iL > 0 ){
      for(j in 1:iL){
        mS0 <- cbind(mShat[-c(1:j),])
        mSj <- cbind(mShat[-c(c(iT-j+1):iT),])
        mGhatj1 <- crossprod(mS0,mSj)/iT #divide by (iT-j) instead?
        mGhatj2 <- crossprod(mSj,mS0)/iT #divide by (iT-j) instead?
        mIhat <- mIhat + vW[j] * (mGhatj1 + mGhatj2)
      }
    }

    ##vcov:
    iN <- length(object$residuals) #divide by n for finite sample version
    object$vcov <- (Jinv %*% mIhat %*% Jinv)/iN
    comment(object$vcov) <- "hac"
    
  } #close if("hac")

  ##return
  ##------
  return(object$vcov) 

} #close vcov.garchx


####################################################
## 4 ADDITIONAL FUNCTIONS
####################################################

##==================================================
## differentiate variable:
gdiff <- function(x, lag=1, pad=TRUE, pad.value=NA)
{
  #check arguments:
  if(lag < 1) stop("lag argument cannot be less than 1")

  #zoo-related:
  iszoo.chk <- is.zoo(x)
  x <- as.zoo(x)
  x.index <- index(x)
  x <- coredata(x)
  isvector <- is.vector(x)
  x <- cbind(x)
  x.n <- NROW(x)
  x.ncol <- NCOL(x)

  #do the differencing:
  xdiff <- x - glag(x, k=lag, pad=TRUE, pad.value=NA)

  #pad operations:
  if(!pad){
    xdiff <- na.trim(as.zoo(xdiff))
    xdiff <- coredata(xdiff)
  }else{
    whereisna <- is.na(xdiff)
    xdiff[whereisna] <- pad.value
  }

  #transform to vector?:
  if(x.ncol==1 && isvector==TRUE){
    xdiff <- as.vector(xdiff)
  }

  #if(is.zoo(x)):
  if(iszoo.chk){
    if(pad){
      xdiff <- zoo(xdiff, order.by=x.index)
    }else{
      xdiff <- zoo(xdiff, order.by=x.index[c(lag+1):x.n])
    } #end if(pad)
  } #end if(iszoo.chk)

  #out:
  return(xdiff)
} #end gdiff

##==================================================
## lag variable k times:
glag <- function(x, k=1, pad=TRUE, pad.value=NA)
{
  #check arguments:
  if(k < 1) stop("Lag order k cannot be less than 1")

  #zoo-related:
  iszoo.chk <- is.zoo(x)
  x <- as.zoo(x)
  x.index <- index(x)
  x <- coredata(x)
  isvector <- is.vector(x)
  x <- cbind(x)
  x.n <- NROW(x)
  x.ncol <- NCOL(x)

  #do the lagging:
  x.nmink <- x.n - k
  xlagged <- matrix(x[1:x.nmink,], x.nmink, x.ncol)
  if(pad){
    xlagged <- rbind( matrix(pad.value,k,x.ncol) , xlagged)
  }

  #transform to vector?:
  if(x.ncol==1 && isvector==TRUE){
    xlagged <- as.vector(xlagged)
  }

  #if(is.zoo(x)):
  if(iszoo.chk){
    if(pad){
      xlagged <- zoo(xlagged, order.by=x.index)
    }else{
      xlagged <- zoo(xlagged, order.by=x.index[c(k+1):x.n])
    } #end if(pad)
  } #end if(iszoo.chk)

  #out:
  return(xlagged)
} #end glag

##==================================================
## simulate random vectors from multivariate normal;
## a speed-optimised version of the rmnorm function
## from the mnormt package of Azzalini (2013),
## http://CRAN.R-project.org/package=mnormt
rmnorm <- function (n, mean = NULL, vcov = 1) 
{
  d <- NCOL(vcov)
  y <- matrix(rnorm(n * d), d, n)
  y <- crossprod(y, chol(vcov))
  if (!is.null(mean)) {
    y <- t(mean + t(y))
  }
  return(y)
}


##==================================================
## perform t-tests under nullity:
ttest0 <- function(x, k=NULL)
{
  ##based on Francq and Thieu (2019), pp. 49-50

  ##check if class is correct:
  if( !is(x, "garchx") ){ stop("'x' not of class 'garchx'") }
  
  ##prepare:
  coefs <- coef.garchx(x)
  coefsNames <- names(coefs)
  n <- length(coefs)
  if(is.null(k)){ k <- 2:n }
  mSigmahat <- vcov.garchx(x) #in fact, mSigmahat/T
  
  ##make matrix w/tests:
  result <- matrix(NA, length(k), 4)
  colnames(result) <- c("coef", "std.error", "t-stat", "p-value")
  rownames(result) <- coefsNames[k]
  result[,"coef"] <- coefs[k]
  
  ##fill matrix w/tests:
  for(i in 1:NROW(result)){
  
    evector <- rep(0, n)
    evector[ k[i] ] <- 1
    stderror <-
      as.numeric(sqrt( rbind(evector) %*% mSigmahat %*% cbind(evector) ))
    statistic <- coefs[k[i]]/stderror
    #note: multiplication by sqrt(T) not needed,
    #since mSigmahat = mSigmahat/n
    
    result[i,"std.error"] <- stderror
    result[i,"t-stat"] <- statistic
    result[i,"p-value"] <- pnorm(statistic, lower.tail=FALSE)
  
  }

  ##result:
  return(result)

} #close ttest0


##==================================================
## perform waldtest under nullity:
waldtest0 <- function(x, r=0, R=NULL, level=c(0.1,0.05,0.01),
  vcov.type=NULL, quantile.type=7, n=20000)
{
  if( !is(x, "garchx") ){ stop("'x' not of class 'garchx'") }
##OLD:
##  if(class(x)!="garchx"){ stop("'x' not of class 'garchx'") }
  
  ##coefs:
  coefs <- as.numeric(coef.garchx(x))
  nmin1 <- length(coefs)-1
  
  ##combination matrix R:
  if(is.null(R)){
    R <- matrix(0, nmin1, nmin1)
    diag(R) <- 1
    R <- cbind( rep(0,nmin1), R)
  }else{
    R <- as.matrix(R)
  }
  rankR <- qr(R)$rank
  if( rankR != NROW(R)){ stop("'R' does not have full rank") }

  ##modify r?
  if( NROW(R)>length(r) ){ r <- rep(r[1], NROW(R)) }
  r <- cbind(r)
  
  ##statistic
  ##---------
  
  SigmahatT <- vcov.garchx(x, vcov.type=vcov.type)
  term3 <- R %*% coefs - r
  term1 <- t(term3)
  term2 <- solve(R %*% SigmahatT %*% t(R))
  statistic <- as.numeric(term1 %*% term2 %*% term3)
  
  ##critical values
  ##---------------
  
  Sigmahat <- nobs(x)*SigmahatT
  Zhat <- t(rmnorm(n, vcov=Sigmahat))
  RZhat <- matrix(NA, NROW(R), n)
  normRZhat <- rep(NA, n) #norm of RZhat
  for(i in 1:n){
    RZhat[,i] <- R %*% Zhat[,i]
    normRZhat[i] <- sum(RZhat[,i]^2)
  }
  critical.values <-
    quantile(normRZhat, prob=1-level, type=quantile.type)
  names(critical.values) <- paste0(100*level,"%")
    
  ##result
  ##------
  
  result <- list(statistic=statistic, critical.values=critical.values)
  return(result)

} #close waldtest0
