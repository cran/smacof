## exploring initial configurations
icExplore <- function(delta, nrep = 100, returnfit = FALSE, ndim = 2, type = c("ratio", "interval", "ordinal","mspline"), 
                       weightmat = NULL, ties = "primary",	verbose = FALSE, 
                       relax = FALSE, modulus = 1, itmax = 1000, eps = 1e-6, spline.degree = 2, spline.intKnots = 2)
{
  ## sanity checks
  type <- match.arg(type, c("ratio", "interval", "ordinal","mspline"), several.ok = FALSE)
  diss <- delta
  if ((is.matrix(diss)) || (is.data.frame(diss))) {
    diss <- strucprep(diss)  #if data are provided as dissimilarity matrix
    attr(diss, "Labels") <- rownames(delta)
  }
  checkdiss(diss)   
  
  n <- attr(diss, 'Size') 
  
  v.stress <- vector() 
  configs <- list()
  simi <- matrix(0, nrow = nrep, ncol = nrep)
  labels <- as.character(1:nrep)
  restot <- list()
  
  for (i in 1:nrep) { 
    if (verbose) cat("IC:", i, "\n")
    z <- matrix(runif(n*ndim, min=0, max=1), nrow=n, ncol=ndim)  ## random inits
    aus1 <- mds(diss, type = type, ndim = ndim, init = z, itmax = itmax, weightmat = weightmat, ties = ties,	
                verbose = FALSE,  relax = relax, modulus = modulus, eps = eps, 
                spline.degree = spline.degree, spline.intKnots = spline.intKnots)
    configs[[i]] <- aus1$conf
    v.stress[i] <- aus1$stress
    restot[[i]] <- aus1 
  }
  for (i in 2:nrep){ 
    je <- i-1
    for (j in 1:je) { 
      a <- configs[[i]]
      b <- configs[[j]]
      aus <- Procrustes(a, b)
      a1 <- c(a)
      b1 <- c(aus$Yhat)
      simi[i,j] <- cor(a1,b1) 
    } 
  }
  sim2 <- simi + t(simi) + diag(nrep) 
  diss1 <- sim2diss(sim2, method = "corr")
  aus2 <- mds(diss1, type = "interval", itmax = itmax)
  x <- aus2$conf[,1]
  y <- aus2$conf[,2]
  
  if (!returnfit) restot <- NULL
  
  res <- list(mdsfit = restot, stressvec = v.stress, conf = cbind(x, y), nrep = nrep, call = match.call())
  class(res) <- "icexplore"
  return(res)
}
