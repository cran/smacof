print.smacofPerm <- function(x,...)
{
  cat("\nCall: ")
  print(x$call)
  cat("\n")
  cat("SMACOF Permutation Test\n")
  cat("Number of objects:",x$nobj,"\n")
  cat("Number of replications:",x$nrep,"\n")
  cat("\n")
  cat("Nonmetric stress value:",x$stress.nm,"\n")
  cat("P-value:",x$pval,"\n")
  cat("\n")
}

