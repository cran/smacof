#Sphere projection, dual algorithm

smacofSphere.dual <- function(delta, penalty = 100, ndim = 2, weightmat = NULL,
                              init = NULL, metric = TRUE, ties = "primary", verbose = FALSE,
                              relax = 1, modulus = 1, itmax = 1000, eps = 1e-6)
{
# penalty ... penalty term kappa >0, 100 is reasonable

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

 mn <- c(1,rep(0,n))
 diss <- as.dist(rbind(0,cbind(0,as.matrix(diss))))         #add row/column

 wgths1 <- as.dist(rbind(0,cbind(0,as.matrix(wgths))))      #distances weights
 wgths2 <- as.dist(outer(mn,mn,function(x,y) abs(x-y)))    
 dhat1 <- as.dist(rbind(0,cbind(0,as.matrix(dhat))))        #0's in the first column
 dhat2 <- mean(sqrt(rowSums(x^2)))*wgths2

 x <- rbind(0,x)
 w <- vmat(wgths1+penalty*wgths2); v<-myGenInv(w); itel<-1;
 d <- dist(x)
 lb <- sum(wgths1*d*dhat1)/sum(wgths1*d^2)
 x <- lb*x
 d <- lb*d
 sold1 <- sum(wgths1*(dhat1-d)^2)
 sold2 <- sum(wgths2*(dhat2-d)^2)
 sold <- sold1+penalty*sold2

 #---------------- begin majorization ---------------- 
 repeat
 {
	b <- bmat(dhat1,wgths1,d)+penalty*bmat(dhat2,wgths2,d)
        y <- v%*%b%*%x
	y <- x+relax*(y-x)
        e <- dist(y)
	ssma1 <- sum(wgths1*(dhat1-e)^2)                       #stress for delta
        ssma2 <- sum(wgths2*(dhat2-e)^2)                       #penalty term 
	ssma <- ssma1+penalty*ssma2                            #joint stress value

        #---- nonmetric MDS --------
        if ((!metric) && ((itel%%modulus) == 0)) {
		if (ties=="primary") daux<-monregP(diss,e,wgths1)
		if (ties=="secondary") daux<-monregS(diss,e,wgths1)
		if (ties=="tertiary") daux<-monregT(diss,e,wgths1)
		daux<-vecAsDist(daux); dhat1<-normDiss(daux,wgths1)
		}
     
	dhat2 <- mean(e[1:n])*wgths2
	snon1 <- sum(wgths1*(dhat1-e)^2)
        snon2 <- sum(wgths2*(dhat2-e)^2)                       
	snon <- snon1+penalty*snon2                            #nonmetric joint stress

        if (verbose) cat("Iteration: ",formatC(itel,width=3, format="d"),"\n",
		" StressOld: ",formatC(c(sold1,sold2,sold),digits=8,width=12,format="f"),"\n",
		" StressSma: ",formatC(c(ssma1,ssma2,ssma),digits=8,width=12,format="f"),"\n",
		" StressNon: ",formatC(c(snon1,snon2,snon),digits=8,width=12,format="f"),"\n",
		"\n\n")

        if (((sold-snon)<eps) || (itel == itmax)) break()      #convergence
        
	x <- y                               #updates
        d <- e
        sold <- snon
        itel <- itel+1
   }
   #-------------- end majorization ---------------

colnames(y) <- paste("D",1:(dim(y)[2]),sep="")
rownames(y) <- labels(diss)
attr(dhat1, "Labels") <- labels(diss)
attr(e, "Labels") <- labels(diss)

if (metric) snon <- NULL          #no non-metric stress
if (!metric) ssma <- NULL

 
result <- list(obsdiss1 = dhat1, obsdiss2 = dhat2, confdiss = e, conf = y, stress.m = ssma, stress.nm = snon,
               ndim = p, model = "Spherical SMACOF (dual)", niter = itel, nobj = n, metric = metric, call = match.call())
class(result) <- c("smacofSP", "smacof")
result
}
