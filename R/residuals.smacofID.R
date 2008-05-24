#residuals of distances for smacof objects

residuals.smacofID <- function(object, ...)
{
  reslist <- list(NULL)
  for (i in 1:length(object$obsdiss)) reslist[[i]] <- (as.matrix(object$obsdiss[[i]] - object$confdiss[[i]]))
  names(reslist) <- names(object$obsdiss)
  return(reslist)  
}
