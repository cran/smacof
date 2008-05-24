`summary.smacofID` <-
function(object, ...)
{
  cat("\n")
  cat("Group Stimulus Space (Joint Configurations):\n")
  print(round(object$gspace, 4))
  #print(round(object$conf,4))
  #cat("\nConfiguration dissimilarities: \n")
  #print(lapply(object$confdiss, round, 4))
}


