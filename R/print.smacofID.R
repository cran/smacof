`print.smacofID` <-
function(x,...)
{
  cat("\n")
  cat("Model:",x$model,"\n")
  cat("Number of objects:",x$nobj,"\n")
  if (!is.null(x$stress.nm)) cat("\nNonmetric stress:",x$stress.nm,"\n")
  if (!is.null(x$stress.m)) cat("\nMetric stress:",x$stress.m,"\n")
  if (!is.null(x$stress.c)) cat("\nConstrained stress",x$stress.co,"\n")
  cat("\nUnconstrained stress",x$stress.uc,"\n")
  cat("Number of iterations:",x$niter,"\n")
  cat("\n")
}


