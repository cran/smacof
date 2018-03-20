`print.smacofR` <-
function(x,...)
{
  cat("\nCall: ")
  print(x$call)
  cat("\n")
  cat("Model:              ",x$model,"\n")
  cat("Number of subjects: ",x$nind,"\n")
  cat("Number of objects:  ",x$nobj,"\n")
  cat("Transformation:     ",x$trans, "\n")
  cat("Conditionality:     ",x$conditionality, "\n")
  
  cat("\nStress-1 value:   ",round(x$stress, 6),"\n")
  cat("Penalized Stress: ", round(x$pstress, 6),"\n")
  cat("Number of iterations:",x$niter,"\n")
  cat("\n")
}

