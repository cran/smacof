`summary.smacofID` <-
function(object, ...)
{
  cat("\n")
  cat("Group Stimulus Space (Joint Configurations):\n")
  print(round(object$gspace, 4))

  cat("\n\n")
  cat("Stress per point:\n")

  spp.perc <- object$spp/sum(object$spp)*100
  sppmat <- cbind(sort(object$spp), sort(spp.perc))
  colnames(sppmat) <- c("SPP","SPP(%)")
  print(round(sppmat, 4))
  cat("\n")

}


