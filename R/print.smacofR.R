`print.smacofR` <-
function(x,...)
{
  cat("\n")
  cat("Model:",x$model,"\n")
  cat("Number of subjects:",x$nind,"\n")
  cat("Number of objects:",x$nobj,"\n")
  
  cat("\nFinal stress value:",x$stress,"\n")
  cat("Number of iterations:",x$niter,"\n")
  cat("\n")
}

