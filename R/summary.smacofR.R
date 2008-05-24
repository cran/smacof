`summary.smacofR` <-
function(object, ...)
{
  cat("\n")
  cat("Subjects configurations (rows):\n")
  print(round(object$conf.row,4))
  cat("\n")
  cat("Objects configurations (columns):\n")
  print(round(object$conf.col,4))
  #cat("\nConfiguration dissimilarities: \n")
  #print(round(object$confdiss,4))
}



