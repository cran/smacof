print.imds <- function(x, ...) {
  cat("\nCall: ")
  print(x$call)
  cat("\n")
  cat("Number of dissimilarity matrices found:", length(x$dissmat),"\n\n")
}