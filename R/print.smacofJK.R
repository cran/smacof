print.smacofJK <- function(x,...)
{
  cat("\nCall: ")
  print(x$call)
  cat("\n")
  cat("SMACOF Jackknife\n")
  cat("Number of objects:",x$nobj,"\n")
  cat("Value loss function:", round(x$loss, 4), "\n")
  cat("Number of iterations:",x$niter,"\n")
  cat("\n")
  cat("Stability measure:",round(x$stab, 4),"\n")
  cat("Cross validity:",round(x$cross,4),"\n")
  cat("Dispersion:",round(x$disp,4),"\n")
  cat("\n")
}

