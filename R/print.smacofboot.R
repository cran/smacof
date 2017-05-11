print.smacofboot <- function(x,...)
{
  cat("\nCall: ")
  print(x$call)
  cat("\n")
  cat("SMACOF Bootstrap:\n")
  cat("Number of objects:",x$nobj,"\n")
  cat("Number of replications:",x$nrep,"\n")
  cat("\n")
  cat("Mean bootstrap stress: ", round(mean(x$stressvec), 4), "\n")
  cat("Stress percentile CI:\n")
  print(round(x$bootci, 4))
  cat("\n")
}

