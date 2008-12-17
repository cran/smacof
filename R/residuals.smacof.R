#residuals of distances for smacof objects

residuals.smacof <- function(object, ...)
{
  # object of class smacof
  return(as.matrix(object$obsdiss - object$confdiss))
}
