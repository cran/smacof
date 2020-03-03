biplotmds.smacof <- function(object, extvar, scale = TRUE) {
  
  if (any(class(object) == "smacofID")) X <- object$gspace else X <- object$conf
  p <- ncol(X)
  
  #if (is.data.frame(extvar)) 
  extvar <- as.data.frame(extvar)
  ## sanity check
  if (nrow(extvar) != nrow(X)) step("Number of rows in extvar needs to match number of objects in configuration!")
  
  rownames(extvar) <- rownames(object$conf)
  ext <- scale(extvar, scale = scale)
  
  regfit <- lm(ext ~ -1 + X)
  regfit$coefficients <- as.matrix(regfit$coefficients)
  if (is.null(colnames(regfit$coefficients))) colnames(regfit$coefficients) <- colnames(extvar)
  rownames(regfit$coefficients) <- colnames(X)
  regsum <- summary(regfit)
  if (ncol(ext) == 1) R2vec <- regsum$r.squared else R2vec <- sapply(regsum, `[[`, "r.squared")
  #names(R2vec) <- gsub("Response ", "", names(R2vec))
  names(R2vec) <- colnames(ext)
  regfit$R2vec <- R2vec
  class(regfit) <- c("mdsbi", "mlm", "lm")
  return(regfit)
}

biplotmds.smacofID <- biplotmds.smacof 