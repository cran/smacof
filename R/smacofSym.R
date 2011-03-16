`smacofSym` <-
function(delta, ndim = 2, weightmat = NULL, init = NULL, 
                      metric = TRUE, ties = "primary",	verbose = FALSE, 
                      relax = FALSE, modulus = 1, itmax = 1000, eps = 1e-6)  
{
# delta ... dissimilarity matrix 
# wghts ... weight structure. if not specified, weights is 1-structure
# p ... number of dimensions
# init ... matrix with starting values of dimension n \times p
# metric ... if TRUE, metric MDS, if FALSE, non-metric
# ties ... ties for pava (primary, secondary, tertiary)
# relax ... relaxation factor
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
  if (p > (n - 1)) stop("Maximum number of dimensions is n-1!")

  nn <- n*(n-1)/2
  m <- length(diss)
  
  if (is.null(attr(diss, "Labels"))) attr(diss, "Labels") <- paste(1:n)
  
  if (is.null(weightmat)) {
    wgths <- initWeights(diss)
  }  else  wgths <- as.dist(weightmat)
  
  #dhat <- normDiss(diss,wgths)            #normalize dissimilarities
  dhat <- normDissN(diss, wgths, 1)        #normalize to n(n-1)/2
                                 
  if (is.null(init)) x <- torgerson(sqrt(diss), p=p) else x <- as.matrix(init)   # x as matrix with starting values   
  if (relax) relax <- 2 else relax <- 1 
  
  w <- vmat(wgths)                        #matrix V of weights and unit vectors
  v <- myGenInv(w)                        #Moore-Penrose inverse
  itel <- 1                               #iteration number
  d <- dist(x)                            #Euclidean distances d(X)
  lb <- sum(wgths*d*dhat)/sum(wgths*d^2)  #denominator: normalization tr(X'VX); 
  x <- lb*x                               #modify x with lb-factor
  d <- lb*d                               #modify d with lb-factor

  sold <- sum(wgths*(dhat-d)^2)           #stress (to be minimized in repeat loop)

 #--------------- begin majorization --------------------
 repeat {                                #majorization loop             
	b <- bmat(dhat,wgths,d)            
        y <- v%*%b%*%x                    #apply Guttman transform denoted as \bar(Y) in the paper
	y <- x+relax*(y-x)                #n \times p matrix of Guttman transformed distances x's
        e <- dist(y)                      #new distance matrix for Y
	ssma <- sum(wgths*(dhat-e)^2)     #stress metric

	if (!metric) {                    #for non-metric MDS only (PAVA)
	 if ((itel%%modulus) == 0) {   
		  if (ties=="primary") daux <- monregP(diss,e,wgths)        #PAVA stuff
		  if (ties=="secondary") daux <- monregS(diss,e,wgths)
		  if (ties=="tertiary") daux <- monregT(diss,e,wgths)
		  dhat <- normDissN(daux,wgths,1)
   }
  }

  snon <- sum(wgths*(dhat-e)^2)     #stress non-metric

  #print out intermediate stress
  if (verbose) cat("Iteration: ",formatC(itel,width=3, format="d")," Stress (not normalized):",
		formatC(c(snon),digits=8,width=12,format="f"),"\n")

  if (((sold-snon)<eps) || (itel == itmax)) break()
  x <- y                           #update configurations
  d <- e                           #update configuration distances
  sold <- snon                     #update stress
  itel <- itel+1	                 #increase iterations
 }
 #------------------ end majorization --------------- 
 
 snon <- snon/nn                   #stress normalization
 ssma <- ssma/nn
 
 if (metric) snon <- NULL          #no non-metric stress
 if (!metric) ssma <- NULL

 colnames(y) <- paste("D",1:(dim(y)[2]),sep="")
 rownames(y) <- labels(diss)
 dhat <- structure(dhat, Size = n, call = quote(as.dist.default(m=b)), class = "dist", Diag = FALSE, Upper = FALSE) 
 attr(dhat, "Labels") <- labels(diss)
 attr(e, "Labels") <- labels(diss)

 confdiss <- normDissN(e, wgths, 1)        #final normalization to n(n-1)/2

 # point stress 
 resmat <- as.matrix(dhat - confdiss)^2    #point stress
 spp <- colMeans(resmat)

 if (itel == itmax) warning("Iteration limit reached! Increase itmax argument!") 
  
#return configurations, configuration distances, normalized observed distances 
result <- list(delta = diss, obsdiss = dhat, confdiss = confdiss, conf = y, stress.m = ssma, stress.nm = snon, spp = spp,
               ndim = p, model = "Symmetric SMACOF", niter = itel, nobj = n, metric = metric, call = match.call()) 
class(result) <- c("smacofB","smacof")
result 
}

