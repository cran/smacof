print.vmu <- function(x,...) {
  cat("\nCall: ")
  print(x$call)
  cat("\n")
  cat("Number of subjects: ", nrow(x$conf.row),"\n", sep = "")
  cat("Number of objects: ", nrow(x$conf.col),"\n", sep = "")
  cat("Number of dimensions: ", x$ndim, "\n", sep = "")
  cat("Variance accounted for: ", round(x$VAF*100, 2), "%\n", sep = "")
  cat("\n")
}
