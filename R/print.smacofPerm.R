print.smacofPerm <- function(x,...)
{
  cat("\nCall: ")
  print(x$call)
  cat("\n")
  cat("SMACOF Permutation Test\n")
  cat("Number of objects:",x$nobj,"\n")
  cat("Number of replications (permutations):",x$nrep,"\n")
  cat("\n")
  cat("Observed stress value:",round(x$stress.obs, 3),"\n")
  if (x$pval == 0) cat("p-value: <0.001 \n") else cat("p-value:", x$pval,"\n")
  cat("\n")
}

