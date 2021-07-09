refit.garchx <- function(object, newy=NULL, newxreg=NULL,
  backcast.values=NULL, reestimate=FALSE, ...)
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
      xreg=newxreg, backcast.values=backcast.values, ...)
  }else{
    coefs <- coef.garchx(object)
    result <- garchx(newy, arch=archArg, garch=garchArg, asym=asymArg,
      xreg=newxreg, initial.values=coefs, backcast.values=backcast.values,
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