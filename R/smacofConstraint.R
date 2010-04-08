#smacof with linear constraints on the configuration (de Leeuw & Heiser, 1980; Borg & Groenen, p. 236)

smacofConstraint <- function(delta, constraint = "linear", external, ndim = 2, weightmat = NULL, startconf = NULL, metric = TRUE,
                             ties = "primary", verbose = FALSE, modulus = 1, itmax = 1000, eps = 1e-6)
{
# diss ... dissimilarity matrix
# constraint ... either "linear", "unique", "diagonal", or a user-specified function
# external ... external data for X-decomposition (Z in paper), or list with "simplex", or "circumplex"
# weightmat ... weight structure. if not specified, weights is 1-structure
# ndim ... number of dimensions
# startconf ... starting configuration
# metric ... if TRUE, metric MDS, if FALSE, non-metric
# ties ... ties for pava (primary, secondary, tertiary)
# modulus ... modulus for nonmetric update
# itmax ... maximum number of iterations
# eps ... change in loss function

  diss <- delta
  if ((is.matrix(diss)) || (is.data.frame(diss))) {
    diss <- strucprep(diss)  #if data are provided as dissimilarity matrix
    attr(diss, "Labels") <- rownames(delta)
  }
  p <- ndim                                     
  n <- attr(diss,"Size")
  nn <- n*(n-1)/2
  m <- length(diss)
  
  if (is.null(attr(diss, "Labels"))) attr(diss, "Labels") <- paste(1:n)
 
  #---- external specification -----
  if (is.data.frame(external)) {
    external <- as.matrix(external)
  }

  if (is.list(external)) {                     
    if (external[[1]] == "simplex") {                           #simplex specification
      d2 <- external[[2]]
      if (d2 >= n) stop("Simplex dimension must be < n!")
      external <- diag(1, n)[,1:d2]
      external[lower.tri(external)] <- 1
    }
    if (external[[1]] == "circumplex") {                        #circumplex specification
      d2 <- external[[2]]
      if (d2 >= n) stop("Circumplex dimension must be <= n!")
      k1 <- external[[3]]
      k2 <- external[[4]]
      if (k2 <= k1) stop("k2 must be > k1")
      external <- matrix(0, nrow = n, ncol = d2)
      ind.all <- expand.grid(1:n, 1:d2)
      inddiff <- apply(ind.all, 1, function(xx) abs(diff(xx)))      
      ind.good <- which((inddiff >= k1) + (inddiff <= k2) == 2) #k1 <= |i-s| <= k2
      el1 <- as.matrix(ind.all[ind.good,])
      external[el1] <- 1
    }
  }

  K <- dim(external)[2]  
  #-------- end external -----------
  
  if (is.null(weightmat)) {
    wgths <- initWeights(diss)
  }  else  wgths <- weightmat
  
 
  dhat <- normDissN(diss,wgths,1)                  #normalize dissimilarities
  w <- vmat(wgths)                              #matrix V
  v <- myGenInv(w)                              #Moore-Penrose inverse
  itel <- 1
  xstart <- startconf

  #----------- pre-specified functions for constraints -----------
  # linear constraint (de Leeuw & Heiser, 1980, p.515), X=ZC 
 if (!is.function(constraint)) { 
  if (constraint == "linear") {
    constrfun <-function(x,w,external) {
      return(external%*%solve(crossprod(external,w%*%external),crossprod(external,w%*%x)))
    }
    if (is.null(xstart)) xstart <- matrix(rnorm(n*p),n,p)                     #starting value for X (before constraints)
  }
  
  # C restricted to be diagonal   
  if (constraint == "diagonal") {
    constrfun <- function(x,w,external) {
      return(external%*%diag(colSums(external*(w%*%x))/colSums(external*(w%*%external))))  
    }
    if (is.null(xstart)) xstart <- matrix(rnorm(n*K),n,K)                     #starting value for X (before constraints) of dimension n x K
  }
  
    
  # X with uniqueness coordinates  
  if (constraint == "unique") {
    constrfun <- function(x,w,external) {
      n <- dim(x)[1]
      p <- dim(x)[2]-n 
      return(cbind(x[,1:p],diag(diag(w%*%x[,p+(1:n)])/diag(w))))
    }
    if (is.null(xstart)) xstart <- cbind(matrix(rnorm(n*p),n,p), diag(1,n))   #starting value for X (before constraints) including diagonal matrix
  }
 } else {   # user-specified
   constrfun <- constraint
   if (is.null(xstart)) stop("Starting configuration must be specified!")
 }
    
  #---------- end pre-specified functions for constraints -------

  
  x <- constrfun(xstart,w,external)                    #compute X

  d <- dist(x)                                         #distances X
  lb <- sum(wgths*d*dhat)/sum(wgths*d^2)               #denominator: normalization tr(X'VX) 
  x <- lb*x                                            #modify x with lb-factor
  d <- lb*d                                            #modify d with lb-factor
  sold <- sum(wgths*(dhat-d)^2)                        #initial stress

  #------- begin majorization -----------
  repeat {                                             #majorization iterations
	b <- bmat(dhat,wgths,d)                        #B matrix
        y <-v%*%b%*%x                                  #Y computation
	y <- constrfun(y,w,external)                   #update Y with corresponding constraints            
        e <- dist(y)                                   #Y distances
	ssma <- sum(wgths*(dhat-e)^2)                  #new stress

	if (!metric) {                                 #nonmetric MDS
	    if ((itel%%modulus) == 0) {
			if (ties=="primary") daux <- monregP(diss,e,wgths)
			if (ties=="secondary") daux <- monregS(diss,e,wgths)
			if (ties=="tertiary") daux <- monregT(diss,e,wgths)
			dhat<-normDissN(daux,wgths,1)
			}
	}
	snon <- sum(wgths*(dhat-e)^2)                  #nonmetric stress
        
	if (verbose) cat("Iteration: ",formatC(itel,width=3, format="d")," Stress: ",
		formatC(c(sold,ssma,snon),digits=8,width=12,format="f"),"\n")

	if (((sold-snon)<eps) || (itel == itmax)) break()   #convergence 

        x <- y                                         #updates
        d <- e
        sold <- snon
        itel <- itel+1	
  }
  #------- end majorization -----------
 
  snon <- snon/nn                   #stress normalization
  ssma <- ssma/nn
  
  if (metric) snon <- NULL          #no non-metric stress
  if (!metric) ssma <- NULL
  
  if (any(is.na(y))) {              #reduce ndim for external == simplex
    csy <- colSums(y)
    ind <- which(is.na(csy))
    y <- y[,-ind]
  }

  colnames(y) <- paste("D",1:(dim(y)[2]),sep="")
  rownames(y) <- labels(diss)
  dhat <- structure(dhat, Size = n, call = quote(as.dist.default(m=b)), class = "dist", Diag = FALSE, Upper = FALSE) 
  attr(dhat, "Labels") <- labels(diss)
  attr(e, "Labels") <- labels(diss)
  
   confdiss <- normDissN(e, wgths, 1)        #final normalization to n(n-1)/2

#return configurations, configuration distances, normalized observed distances 
  result <- list(obsdiss = dhat, confdiss = confdiss, conf = y, stress.m = ssma, stress.nm = snon,
               ndim = p, model = "SMACOF constraint", niter = itel, nobj = n, metric = metric, call = match.call()) 
  class(result) <- c("smacofB","smacof")
  result 
}

