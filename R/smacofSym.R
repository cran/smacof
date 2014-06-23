`smacofSym` <- function(delta, ndim = 2, type = c("ratio", "interval", "ordinal","mspline"), 
                        weightmat = NULL, init = NULL, ties = "primary",	verbose = FALSE, 
                        relax = FALSE, modulus = 1, itmax = 1000, eps = 1e-6, 
                        spline.degree = 2, spline.intKnots = 2)  {
  # delta ... dissimilarity matrix 
  # wghts ... weight structure. if not specified, weights is 1-structure
  # p ... number of dimensions
  # init ... matrix with starting values of dimension n \times p
  # type ... either "ratio", "interval", "ordinal", "spline" (replaces metric)
  # ties ... ties for pava (primary, secondary, tertiary)
  # relax ... relaxation factor
  # modulus ... modulus for nonmetric update
  # itmax ... maximum number of iterations
  # eps ... change in loss function
  # spline.degree ... degree of the spline in case a spline transformation is chosen
  # spline.intKnots ... number of interior knots for the spline in case a spline transformation is chosen
  
  ## --- sanity checks
  type <- match.arg(type, c("ratio", "interval", "ordinal","mspline"), several.ok = FALSE)
  
  diss <- delta
  if ((is.matrix(diss)) || (is.data.frame(diss))) {
    diss <- strucprep(diss)  #if data are provided as dissimilarity matrix
    attr(diss, "Labels") <- rownames(delta)
  }
  
  p <- ndim                                     
  n <- attr(diss,"Size")
  if (p > (n - 1)) stop("Maximum number of dimensions is n-1!")
  
  nn <- n*(n-1)/2
  m <- length(diss)
  
  if (is.null(attr(diss, "Labels"))) attr(diss, "Labels") <- paste(1:n)
  
  if (is.null(weightmat)) {
    wgths <- initWeights(diss)
  }  else  {
    wgths <- as.dist(weightmat)
  }
  
  # Prepare for optimal scaling
  trans <- type
  if (trans=="ratio"){
    trans <- "none"
  } else if (trans=="ordinal" & ties=="primary"){
    trans <- "ordinalp"
  } else if(trans=="ordinal" & ties=="secondary"){
    trans <- "ordinals"
  } else if(trans=="ordinal" & ties=="tertiary"){
    trans <- "ordinalt"
  } else if(trans=="spline"){
    trans <- "mspline"
  }
  disobj <- transPrep(diss,trans = trans, spline.intKnots = spline.intKnots, spline.degree = spline.degree)
  
  #dhat <- normDiss(diss,wgths)            #normalize dissimilarities
  dhat <- normDissN(diss, wgths, 1)        #normalize to n(n-1)/2
  
  if (is.null(init)) x <- torgerson(sqrt(diss), p=p) else x <- as.matrix(init)   # x as matrix with starting values   
  if (relax) relax <- 2 else relax <- 1 
  
  w    <- vmat(wgths)                      #matrix V of weights and unit vectors
  v    <- myGenInv(w)                      #Moore-Penrose inverse
  itel <- 1                                #iteration number
  d    <- dist(x)                          #Euclidean distances d(X)
  lb   <- sum(wgths*d*dhat)/sum(wgths*d^2) #denominator: normalization tr(X'VX); 
  x    <- lb*x                             #modify x with lb-factor
  d    <- lb*d                             #modify d with lb-factor
  
  sold <- sum(wgths*(dhat-d)^2)/nn         #stress (to be minimized in repeat loop)
  
  #--------------- begin majorization --------------------
  repeat {                                #majorization loop             
    b <- bmat(dhat,wgths,d)            
    y <- v%*%b%*%x                    #apply Guttman transform denoted as \bar(Y) in the paper
    y <- x+relax*(y-x)                #n \times p matrix of Guttman transformed distances x's
    e <- dist(y)                      #new distance matrix for Y
    ssma <- sum(wgths*(dhat-e)^2)     #stress metric
    
    dhat2 <- transform(e, disobj, w = wgths, normq = nn)
    dhat <- dhat2$res
    
    snon <- sum(wgths*(dhat-e)^2)/nn     #stress non-metric
    
    #print out intermediate stress
    if (verbose) cat("Iteration: ",formatC(itel,width=3, format="d"),
                     " Stress (normalized):", formatC(c(snon),digits=8,width=12,format="f"),
                     " Difference: ", formatC(sold-snon,digits=8,width=12,format="f"),"\n")
    
    if (((sold-snon)<eps) || (itel == itmax)) break()
    x <- y                           #update configurations
    d <- e                           #update configuration distances
    sold <- snon                     #update stress
    itel <- itel+1	                 #increase iterations
  }
  #------------------ end majorization --------------- 
  
  stress <- sqrt(snon)                   #stress normalization
  
  colnames(y) <- paste("D",1:(dim(y)[2]),sep="")
  rownames(y) <- labels(diss)
  dhat <- structure(dhat, Size = n, call = quote(as.dist.default(m=b)), class = "dist", Diag = FALSE, Upper = FALSE) 
  attr(dhat, "Labels") <- labels(diss)
  attr(e, "Labels") <- labels(diss)
  
  confdiss <- normDissN(e, wgths, 1)        #final normalization to n(n-1)/2
  
  ## point stress 
  resmat <- as.matrix(dhat - confdiss)^2    #point stress
  spp <- colMeans(resmat)
  
  if (itel == itmax) warning("Iteration limit reached! Increase itmax argument!") 
  
  #return configurations, configuration distances, normalized observed distances 
  result <- list(delta = diss, obsdiss = dhat, confdiss = confdiss, iord = dhat2$iord.prim, conf = y, stress = stress, spp = spp,
                 ndim = p, model = "Symmetric SMACOF", niter = itel, nobj = n, type = type, call = match.call()) 
  class(result) <- c("smacofB","smacof")
  result 
}

