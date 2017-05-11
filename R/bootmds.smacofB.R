bootmds.smacofB <- function(object, data,  method.dat = "pearson", nrep = 100, alpha = 0.05, verbose = FALSE, ...) 
{
  ## object... object of class smacofB (from smacofSym, smacofConstraint)
  if (class(object)[1] != "smacofB") stop("Bootstrap is currently implemented for objects of class smacofB from smacofSym() only! \n")
  if (object$model == "SMACOF constraint") stop("Bootstrap is currently implemented for smacofSym() objects only! \n")
  
    
  method.dat <- match.arg(method.dat, c("pearson", "spearman", "kendall", "euclidean", "maximum", "manhattan", "canberra", "binary"))
  n <- object$nobj          ## number of objects
  if (!missing(data)) if(ncol(data) != n) stop("Number of columns need to match number of MDS objects!") 
  p <- object$ndim          ## number of dimensions
  val <- object$stress     
  smacall <- object$call
  
  N <- dim(data)[1]
  coord <- list()
  x <- vector()
  y <- vector()
  stressvec <- c()
  
  for (i in 1:nrep) {
    st <- data[sample(1:N, size = N, replace = TRUE), ]      ## bootstrap sample data
    if (verbose) cat("Replication: ", i, "\n")
    
    ## compute input dissimilarities
    if (any(method.dat == c("pearson", "spearman", "kendall"))) {
      r <- cor(st, method = method.dat, use = "pairwise.complete.obs")    ## compute proximities (correlation)
      r <- sim2diss(r, ...) 
    } else {
      r <- as.matrix(dist(t(st), method = method.dat))                    ## compute dissimilarity
    }
   
    smacall$delta <- r
    o <- eval(smacall)  
    stressvec[i] <- o$stress
    fit <- Procrustes(object$conf, o$conf)
    coord[[i]] <- fit$Yhat
  }  
  
  M <- list()  
  xy <- matrix(NA, nrow = nrep, ncol = object$ndim)
  for (k in 1:n){ 
    for (i in 1:nrep) {
      xy[i,] <- coord[[i]][k,]
    }
    M[[k]] <- cov(xy) 
  } 
  
  names(M) <- attr(object$confdist, "Labels")
  bootci <- quantile(stressvec, probs = c(alpha/2, (1-alpha/2)))
  
	result <- list(cov = M, conf = object$conf, bootconf = coord, stressvec = stressvec, nrep = nrep, nobj = n, alpha = alpha, bootci = bootci, call = match.call())
  class(result) <- "smacofboot"
  result
}
    