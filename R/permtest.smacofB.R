## MDS permutation test
permtest.smacofB <- function(object, nrep = 100, verbose = FALSE)
{
## val ... stress value  
## n... number of objects
## p... number of dimensions
    
    if (class(object)[1] != "smacofB") stop("Permutation test is currenlty implemented for objects of class smacofB from smacofSym() only! \n")
    if (object$model == "SMACOF constraint") stop("Permutation test is currenlty implemented for smacofSym() objects only! \n")
    if (object$metric) stop("Permutation test implemented for nonmetric SMACOF only! \n")
   
    n <- object$nobj          ## number of objects
    p <- object$ndim          ## number of dimensions
    val <- object$stress.nm   ## nonmetric stress
    
    m <- n*(n-1)/2            ## number of lower triangular elements
    irep <- 1
    str <- rep (0, nrep)      ## vector for stress values
    
    repeat 
    {
        delta <- matrix (0, n, n)
        delta[outer(1:n, 1:n, ">")] <- sample (1:m,m)              ## sample dissimilarity matrix
        delta <- delta + t(delta)                                  
        smRes <- smacofSym(delta, ndim = p, metric = FALSE)        ## non-metric smacof
        str[irep] <- smRes$stress.nm                               ## store stress of no-structure matrix 
        
        if (verbose) cat("Replication: ", formatC (irep, digits=3, width=3), "Stress: ", formatC (str[irep], digits=10, width=15, format="f"), "\n")
		
        if (irep == nrep) break()
  
        irep <- irep + 1
    }
    
    pval <- length(which(str < val))/nrep
      
    result <- list(stressvec = str, stress.nm = val, pval = pval, nobj = n, nrep = nrep, call = match.call())
    class(result) <- "smacofPerm"
    result
}