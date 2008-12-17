#Sphere projection, primal algorithm

smacofSphere.primal <- function(delta, ndim = 2, weightmat = NULL, init = NULL,
                              	metric = TRUE, ties = "primary", verbose = FALSE,
                                modulus = 1, itmax = 1000, eps = 1e-6)

{

 diss <- delta
 if ((is.matrix(diss)) || (is.data.frame(diss))) diss <- strucprep(diss)  #if data are provided as dissimilarity matrix
 p <- ndim
 n <- attr(diss,"Size")
 m <- length(diss)
 if (is.null(attr(diss, "Labels"))) attr(diss, "Labels") <- paste(1:n)
  

 if (is.null(weightmat)) {
    wgths <- initWeights(diss)
 }  else  wgths <- weightmat

 dhat <- normDiss(diss,wgths)            #normalize dissimilarities
 if (is.null(init)) x <- torgerson(sqrt(diss), p=p) else x <- init   # x as matrix with starting values


 w <- vmat(wgths)
 v <- myGenInv(w)
 itel<-1;

 x <- x/sqrt(rowSums(x^2))
 #FIXME!!!
 d <- dist(x)                           #distance computation (to be extended with geodesic)

 lb <- sum(wgths*d*dhat)/sum(wgths*d^2)
 x <- lb*x
 d <- lb*d
 sold <- sum(wgths*(dhat-d)^2)          #initial stress

 #------------------------- begin majorization ---------------------------------
 repeat {
 	 b <- bmat(dhat,wgths,d)
   y <- v%*%b%*%x                       #Guttman transform
	 y <- sphereProj(y,w)                 #projection on the sphere
	 #FIXME!!!
   e <- dist(y)                         #distances for Y (to be enhanced with geodesics
	 ssma <- sum(wgths*(dhat-e)^2)        #metric stress

   if (!metric) {                       #nonmetric versions
	    if ((itel%%modulus) == 0) {
			if (ties=="primary") daux <- monregP(diss,e,wgths)
			if (ties=="secondary") daux <- monregS(diss,e,wgths)
			if (ties=="tertiary") daux <- monregT(diss,e,wgths)
			daux <- vecAsDist(daux)
      dhat <- normDiss(daux,wgths)
			}
	}
  snon <- sum(wgths*(dhat-e)^2)        #nonmetric stress
	if (verbose) cat("Iteration: ",formatC(itel,width=3, format="d")," Stress: ",
		formatC(c(sold,ssma,snon),digits=8,width=12,format="f"),"\n")
	if (((sold-snon)<eps) || (itel == itmax)) break()

  x <- y                               #updates
  d <- e
  sold <- snon
  itel <- itel+1
 }
 #----------------------------- end majorization -------------------------------

colnames(y) <- paste("D",1:(dim(y)[2]),sep="")
rownames(y) <- labels(diss)
attr(dhat, "Labels") <- labels(diss)
attr(e, "Labels") <- labels(diss)

if (metric) snon <- NULL          #no non-metric stress
if (!metric) ssma <- NULL


result <- list(obsdiss = dhat, confdiss = e, conf = y, stress.m = ssma, stress.nm = snon,
               ndim = p, model = "Spherical SMACOF (primal)", niter = itel, nobj = n)
class(result) <- c("smacofSP", "smacof")
result
}
