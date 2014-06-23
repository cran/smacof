`print.smacofB` <-
function(x,...)
{
  cat("\nCall: ")
  print(x$call)
  cat("\n")
  cat("Model:",x$model,"\n")
  cat("Number of objects:",x$nobj,"\n")
  cat("Stess-1 value:",x$stress,"\n")
  cat("Number of iterations:",x$niter,"\n")
  cat("\n")
}

