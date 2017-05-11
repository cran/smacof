print.icexplore <- function(x,...)
{
  cat("\nCall: ")
  print(x$call)
  cat("\n")
  cat("Number of replications:",x$nrep,"\n")
  cat("\nBest stress value:", round(min(x$stressvec), 4),"\n")
  cat("Average stress value: ", round(mean(x$stressvec), 4), "\n")
  cat("Stress quantiles:\n")
  print(round(quantile(x$stressvec, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)), 4))
  cat("\n")
}

