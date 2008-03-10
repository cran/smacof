smacofRect <- function(delta, ndim = 2, weightmat = NULL, init = NULL, verbose = FALSE,
                       itmax = 1000, reg = 1e-6, eps = 1e-6)

# init ... either a list of 2 matrices of dimension n \times p and m \times p with 
# starting values. if NULL, svd is used.
{

  diss <- delta
  rnames <- rownames(delta)
  if (is.data.frame(diss)) diss <- as.matrix(diss)
  n <- dim(diss)[1]                       #number of individuals
  m <- dim(diss)[2]                       #number of objects
  p <- ndim

  if (is.null(weightmat)) {
    w <- matrix(1,n,m)                    #initialize weights (as 1)
  } else w <- weightmat
  

  itel <- 1
  delta <- ifelse(is.na(diss),0,diss)     #replace NA's by 0
  #delta <- delta/sqrt(sum(w*delta^2))*sqrt(n*m)       #normalize dissimilarities
  
  delta_plus <- ifelse(delta>=0,delta,0)  #delta decomposition (+)
  delta_min <- ifelse(delta<=0,-delta,0)  #delta decomposition (-) (if all >0 --> complete 0)

  if (is.list(init)) {
    x <-init[[1]]                         #list as input structure
    y <-init[[2]]
  } else {
    e <- delta_plus^2
    e <- -0.5*(e-outer(rowSums(e)/m,colSums(e)/n,"+")+(sum(e)/(n*m)))
    #e <- e/sqrt(sum(e^2))*sqrt(n*m) 
    
    z <- svd(e,nu=p,nv=0)                 #SVD for e (pos. distances)
    x<-z$u                                #starting value for x
    y<-crossprod(e,x)                     #starting value for y
  }

  d <- distRect(x,y,reg)                  #n times m of reproduced diss
  #d <- d/sqrt(sum(w*d^2))*sqrt(n*m)
  lold <- sum(w*(delta-d)^2)              #stress value
  
  #------------------- begin majorization -----------------------------
  repeat {
	ww <- w*(1+(delta_min/d)) 
        wr <- rowSums(ww)
        wc <- colSums(ww)

        v <- solve(diag(wc)+(1/m) - crossprod(ww,ww/wr)) -(1/m)

        b <- w*delta_plus/d                #B matrix
        br <- rowSums(b)                   #rows B
        bc <- colSums(b)                   #columns W

        xraw <- (br*x)-(b%*%y)
        yraw <- (bc*y)-crossprod(b,x)

        y <- v%*%(yraw+crossprod(ww,xraw/wr)) #x update 
        x <- (xraw+(ww%*%y))/wr               #y update

        d <- distRect(x,y,reg)             #compute distances (update)
        lnew <- sum(w*(delta-d)^2)         #compute stress

        if (verbose) {
		cat("Iteration: ",formatC(itel,digits=6,width=6),
	    	"   Stress:   ",formatC(lold,digits=6,width=12,format="f"),
	    	" ==>",formatC(lnew,digits=6,width=12,format="f"),"\n")
	}

        if (((lold-lnew) < eps) || (itel==itmax)) break() 
        
	      lold <- lnew                       #update stress
        itel <- itel+1
  }
  #-------------------- end majorization --------------------------

colnames(y) <- colnames(x) <- paste("D",1:(dim(y)[2]),sep="")
rownames(x) <- rownames(diss) <- rownames(d) <- rnames

 
#return configuration distances, row and column configurations, stress 
result <- list(obsdiss = diss, confdiss = d, conf.row = x, conf.col = y, stress = lnew, 
               ndim = p, model = "Rectangular smacof", niter = itel, nind = n, nobj = m, metric = TRUE, call = match.call()) 
class(result) <- "smacofR"
result 
}
