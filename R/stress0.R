stress0 <- function(delta, init, type = c("interval", "ratio", "ordinal","mspline"), 
                    weightmat = NULL, ties = "primary",	spline.degree = 2, spline.intKnots = 2)  {
  
  ## utility function to determine stress of niter = 0.  
 
  
  ## --- sanity checks
  type <- match.arg(type, c("interval", "ratio", "ordinal","mspline"), several.ok = FALSE)
  ndim <- ncol(init)
  
  diss <- delta
  if ((is.matrix(diss)) || (is.data.frame(diss))) {
    diss <- strucprep(diss)  #if data are provided as dissimilarity matrix
    attr(diss, "Labels") <- rownames(delta)
  }
  checkdiss(diss)           ## check whether dissimilarities are all positive       
  
  p <- ndim                                     
  n <- attr(diss,"Size")
  if (p > (n - 1)) stop("Maximum number of dimensions is n-1!")
  
  nn <- n*(n-1)/2
  m <- length(diss)
  
  if (is.null(attr(diss, "Labels"))) attr(diss, "Labels") <- paste(1:n)
  
  ## --- weights
  if (is.null(weightmat)) {
    wgths <- initWeights(diss)
  }  else  {
    wgths <- as.dist(weightmat)
  }
  
  ## --- starting values
  x <- initConf(init, diss, n, p)
  xstart <- x
  
  ## --- Prepare for optimal scaling
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
  
  ## dhats and missings
  dhat <- normDissN(diss, wgths, 1)        ## normalize to n(n-1)/2
  dhat[is.na(dhat)] <- 1     ## in case of missing dissimilarities, pseudo value for dhat 

  w    <- vmat(wgths)                      #matrix V of weights and unit vectors
  v    <- myGenInv(w)                      #Moore-Penrose inverse
  itel <- 1                                #iteration number
  d    <- dist(x)                          #Euclidean distances d(X)
  lb   <- sum(wgths*d*dhat, na.rm = TRUE)/sum(wgths*d^2) #denominator: normalization tr(X'VX); 
  x    <- lb*x                             #modify x with lb-factor
  d    <- lb*d                             #modify d with lb-factor
  
  #sold <- sum(wgths*(dhat-d)^2, na.rm = TRUE)/nn         #stress (to be minimized in repeat loop)
  
  dhat2 <- transform(d, disobj, w = wgths, normq = nn)  ## dhat update
  dhat <- dhat2$res
  snon <- sum(wgths*(dhat-d)^2)/nn   
  stress <- sqrt(snon)
  return(stress)
}