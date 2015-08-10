permtest.smacofR <- function(object, nrep = 100, verbose = TRUE)
{
## val ... stress value  
## n... number of objects
## p... number of dimensions
## ... additional arguments to be passed from smacofRect  
    
    #if (class(object)[1] != "smacofR") stop("Permutation test is currenlty implemented for objects of class smacofB from smacofSym() only! \n")
    
    data <- object$obsdiss
    m <- object$nobj          ## number of objects (columns)
    n <- object$nind          ## number of observations (rows)
    p <- object$ndim          ## number of dimensions
    val <- object$stress      ## metric stress (normalized)
    smacall <- object$call
    
    stressvec <- rep (0, nrep)      ## vector for stress values
    congmat <- matrix(0, nrep, n)
    #perms <- shuffleSet(m, nset = nper)
    
    for (irep in 1:nrep) {
      permmat <- t(apply(data, 1, function(pp) {     ## computes permuted matrix
        ind <- sample(1:m, m)
        pp[ind]
      }))
      
      smacall$delta <- permmat
      resperm <- eval(smacall) 
      #resperm <- smacofRect(permmat, ndim = p, ...)
      stressvec[irep] <- resperm$stress
      congmat[irep, ] <- resperm$congvec
              
      if (verbose) cat("Permutation: ", formatC (irep, digits=3, width=3), "Stress: ", formatC (stressvec[irep], digits=10, width=15, format="f"), "\n")
    }
    pval <- length(which(stressvec < val))/nrep
      
    result <- list(stressvec = stressvec, stress.obs = val, pval = pval, 
                   congmat = congmat, nobj = n, nrep = nrep, call = match.call())
    class(result) <- "smacofPerm"
    result
}