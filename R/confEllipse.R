confEllipse <- function(object)UseMethod("confEllipse")

confEllipse.smacofID <- function(object) {
  
  ## sanity checks
  if (object$type != "ratio") stop("Confidence ellipses implemented for ratio transformation only!\n")
  
  ## data preparation
  if (any(class(object) == "smacofID")) X <- object$gspace else X <- object$conf      ## configuration
  delta <- normDeltaOff(object$delta, object$weightmat)         ## normalize delta
  w <- object$weightmat                                         ## getting weight matrix in shape
  if (!is.list(w)) w <- list(w)
  wlist <- lapply(w, as.matrix)
  
  
  ## call additional optimization steps
  if (is.null(object$constraint)) object$constraint <- "identity"             ## for ordinary smacof solutions
  if (object$constraint == "identity") { 
    h <- smacofNODIFF(wlist, delta, X)
    h$b <- repList(diag(ncol(X)), length(wlist))                     ## create individual weight matrix
  }  
  if (object$constraint == "diagonal") h <- smacofINDSCAL(wlist, delta, object$cweights, X)
  if (object$constraint == "idioscal") h <- smacofIDIOSCAL(wlist, delta, object$cweights, X)
  
  ## compute derivatives
  hh <- smacofDerivativesX(wlist, delta, h$b, h$x)
  
  rownames(h$x) <- rownames(X)
  colnames(h$x) <- colnames(X)
  result <- list(X = h$x, h = hh$h, s = hh$s)
  class(result) <- "confell"
  return(result)
}

confEllipse.smacofB <- confEllipse.smacofID
