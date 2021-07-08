refit.garchx <- function(object, newy=NULL, newxreg=NULL,
  reestimate=FALSE, ...)
{
  ##check:
  if( is.null(newy)){ stop("'newy' cannot be NULL") }
  if( !is.null(object$xreg) && is.null(newxreg) ){
    stop("'newxreg' is missing")
  }

  ##splice old and new data:
  y <- c(object$y.coredata, coredata(as.zoo(newy)))
  if( is.null(newxreg) ){
    xregArg <- NULL
  }else{
    xregArg <- rbind(object$xreg, newxreg)
  }
  
  ##obtain spec:    
  archArg <- object$arch
  garchArg <- object$garch 
  asymArg <- object$asym

  ##refit:
  if(reestimate){
    result <- garchx(y, arch=archArg, garch=garchArg, asym=asymArg,
      xreg=xregArg)
  }else{
    coefs <- coef.garchx(object)
    result <- garchx(y, arch=archArg, garch=garchArg, asym=asymArg,
      xreg=xregArg, initial.values=coefs, estimate=FALSE, turbo=TRUE)
    result$fitted <- fitted.garchx(result)
    result$residuals <- residuals.garchx(result)
    result$hessian <- object$hessian
    result$vcov <- object$vcov
  }
    
  ##return result:
  return(result)
}