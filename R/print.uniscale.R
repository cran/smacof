print.uniscale <- function(x,...)
  {
    cat("\nCall: ")
    print(x$call)
    cat("\nStress-1 value:", round(x$stress, 4),"\n")
    cat("Number of accepted permutations:",x$npermscale,"\n")
    cat("Number of possible permutations:",x$npermtot,"\n")
    cat("Number of objects:",x$nobj,"\n")
    cat("\n")
  }
