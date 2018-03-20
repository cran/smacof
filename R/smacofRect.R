smacofRect <- function(delta, ndim = 2, type = c("ratio", "interval", "ordinal", "mspline"),
                       conditionality = c("matrix", "row"), lambda = 0.5, omega = 0.1, 
                       circle = c("none", "row", "column"), weightmat = NULL, init = NULL, 
                       ties = c("primary", "secondary"), verbose = FALSE, relax = TRUE, itmax = 10000,  
                       eps = 1e-6, spline.degree = 2, spline.intKnots = 2)
  
  # init ... either a list of 2 matrices of dimension n \times p and m \times p with 
  # starting values. if NULL, svd is used.
{
  circle <- match.arg(circle, c("none", "row", "column"), several.ok = FALSE)
  type <- match.arg(type, c("ratio", "interval", "ordinal", "mspline"), several.ok = FALSE)
  ties <- match.arg(ties, c("primary","secondary"), several.ok = FALSE)
  conditionality <- match.arg(conditionality, c("matrix", "row"), several.ok = FALSE)
  diss <- delta
  rnames <- rownames(delta)
  if (is.data.frame(diss)) diss <- as.matrix(diss)
  checkdiss(diss)

  n <- dim(diss)[1]                       #number of individuals
  m <- dim(diss)[2]                       #number of objects
  p <- ndim

  if (is.null(weightmat)) {
    w <- matrix(1,n,m)                    #initialize weights (as 1)
  } else w <- weightmat
  wpp <- sum(w)

  delta <- ifelse(is.na(diss),0,diss)     #replace NA's by 0

  ## --- Prepare for optimal scaling
  trans <- type
  if (trans=="ratio"){
    trans <- "none"
  } else if (trans=="ordinal" & ties=="primary"){
    trans <- "ordinalp"
  } else if(trans=="ordinal" & ties=="secondary"){
    trans <- "ordinals"
  } else if(trans=="spline"){
    trans <- "mspline"
  }
  disobj <- list()
  if (conditionality == "matrix"){
    disobj[[1]] <- transPrep(as.vector(delta), trans = trans, spline.intKnots = spline.intKnots, spline.degree = spline.degree)
    if (trans == "mspline") disobj[[1]]$base <- cbind(rep(1, nrow(disobj[[1]]$base)), disobj[[1]]$base)
    tt <- transform(as.vector(delta), disobj[[1]], w = as.vector(w)) #, normq = wpp
    dhat <- matrix(tt$res, n, m)  ## dhat update
    tt <- vector(mode = "list", length = 1)
  } else { ## conditionality == "row"
    dhat <- matrix(0, n, m)
    for (i in 1:n) {
      disobj[[i]] <- transPrep(delta[i, ], trans = trans, spline.intKnots = spline.intKnots, spline.degree = spline.degree)
      if (trans == "mspline") disobj[[i]]$base <- cbind(rep(1, nrow(disobj[[i]]$base)), disobj[[i]]$base)
      dhat[i, ] <- transform(delta[i, ], disobj[[i]], w = w[i, ], normq = sum(w))$res  ## dhat update
    }
    tt <- vector(mode = "list", length = n)
  }
  ## --- end optimal scaling prep


  itel <- 1
  #delta <- delta/sqrt(sum(w*delta^2))*sqrt(n*m)       #normalize dissimilarities

  #delta_plus <- ifelse(delta>=0,delta,0)  #delta decomposition (+)
  #delta_min <- ifelse(delta<=0,-delta,0)  #delta decomposition (-) (if all >0 --> complete 0)

  if (is.list(init)) {
    x <-init[[1]]                         #list as input structure
    y <-init[[2]]
  } else {
    e <- dhat^2
    e <- -0.5*(e-outer(rowSums(e)/m,colSums(e)/n,"+")+(sum(e)/(n*m)))
    #e <- e/sqrt(sum(e^2))*sqrt(n*m)

    z <- svd(e,nu=p,nv=0)                 #SVD for e (pos. distances)
    x<-z$u                                #starting value for x
    y<-crossprod(e,x)                     #starting value for y
  }

  if (circle != "none"){
    r <- projCircle(x,y,x,y,circle=circle)
    x <- r$x
    y <- r$y
    wr <- rowSums(w)
    wc <- colSums(w)
    lambda <- 2*max(c(wr,wc))
  }

  d <- distRect(x,y,0)                  #n times m of reproduced diss
  coefOfVar <- function(x, w){            # Compute coefficient of variation
    av <- sum(x * w)/sum(w)
    va <- sum((x - av)^2 * w)/sum(w)
    return(va/av)
  }

  # psychometrika stress
  pstress <- function( dhat, d, w, omega, lambda, wpp ){
    if ( conditionality == "matrix" ){
      c1 <- nrmd2 <- sum( w * ( dhat - d )^2 )
      sumw <- sum( w )
      mnc1 <- nrmd2 / wpp
      wsum <- sum( w * dhat )
      nrmw2 <- sum( w * dhat^2 )
      nrmm2 <- wsum^2 / wpp
      nrmv2 <- nrmw2 - nrmm2
      c2 <- mnc2 <- ( nrmv2 + omega * nrmm2 ) / nrmv2
      g <- sqrt( mnc1^lambda * mnc2 )
    } else { ## conditionality == "row"
      c1 <- nrmd2 <- rowSums( w * ( dhat - d )^2 )
      sumw <- rowSums( w )
      wsum <- rowSums( w * dhat )
      nrmw2 <- rowSums( w * dhat^2 )
      nrmm2 <- wsum^2 / sumw
      nrmv2 <- nrmw2 - nrmm2
      c2 <- ( nrmv2 + omega * nrmm2 ) / nrmv2
      mnc1 <- sum( c1 )
      mnc2 <- sum( c2 )
      g <- sqrt( ( mnc1 / wpp )^lambda * mnc2 )
    }
    return( list( pstress = g, g = g, c1 = c1, c2 = c2, nrmd2 = nrmd2, sumw = sumw, wsum = wsum, nrmw2 = nrmw2, nrmm2 = nrmm2, nrmv2 = nrmv2 ) )
  } # psychometrika pstress

  ps   <- pstress(dhat, d, w, omega, lambda, wpp)   #pstress value
  lold <- ps$pstress

  if (circle == "none") {
    ww <- w
    wr <- rowSums(ww)
    wc <- colSums(ww)
    v <- solve(diag(wc) + (1/m) - crossprod(ww, ww/wr)) - (1/m)
  }


  #------------------- begin majorization -----------------------------
  repeat {
    d.is.0 <- d < eps 
    if (circle == "none") {
      b  <- w * dhat * (!d.is.0) / (d + d.is.0) # B matrix
      br <- rowSums(b)                   #rows B
      bc <- colSums(b)                   #columns W

      xraw <- (br * x) - ( b %*% y)
      yraw <- (bc * y) - crossprod(b, x)


      xold <- x
      yold <- y
      y <- v %*% (yraw + crossprod(ww, xraw/wr)) #x update
      x <- (xraw + (ww %*% y))/wr                #y update
      if (relax & itel > 100){
        x <- 2 * x - xold
        y <- 2 * y - yold
      }
      
    } else {
      b  <- w * (1 - (!d.is.0) * dhat / (d + d.is.0))  #B matrix
      br <- rowSums(b)                   #rows B
      bc <- colSums(b)                   #columns W
      xunc <- x - outer(br, rep(1/lambda, p), "*") * x + b %*% (y/lambda)
      yunc <- y - outer(bc, rep(1/lambda, p), "*") * y + t(b) %*% (x/lambda)
      r <- projCircle(xunc, yunc, x, y, circle = circle)
      x <- r$x
      y <- r$y
    }

    d <- distRect(x,y,0)             #compute distances (update)

    #lnew <- sum(w*(dhat - d)^2)/wpp     #compute stress

    # Update dhats

    dhat.old <- dhat

    # psychometrika target
    if ( conditionality == "matrix" ) {
      nrmw <- sqrt( ps$nrmw2 )
      nrmm <- sqrt( ps$nrmm2 )
      nrmv <- sqrt( ps$nrmv2 )
      g1 <- ps$nrmd2 / wpp
      g2 <- ps$nrmv2 + omega * ps$nrmm2
      g3 <- sqrt( ps$nrmv2 )
      g <- g1^( lambda / 2 ) * sqrt( g2 ) / g3
      alpha2 <- 0.5 * lambda * sqrt( g2 ) * g1^( lambda / 2 - 1 )
      alpha3 <- 0.5 * g1^( lambda / 2 ) / sqrt( g2 )
      fmin <- min( dhat )
      tau1 <- ps$wsum / wpp
      tau2   <- 0.5 / eps
      if ( fmin >= eps ) tau2 <- 0.5 / fmin
      tau4   <- tau1 / nrmv
      beta1  <- 1 / wpp
      beta3  <- 1 + omega + 2 * omega * tau1 * tau2
      beta5  <- tau2 * tau4
      lower  <- alpha2 * beta1 + alpha3 * beta3 + g * beta5
      t1     <- tau2 * as.vector( dhat ) - 0.5
      t2     <- as.vector( dhat ) / ( 2 * nrmv )
      b1     <- as.vector( d ) / wpp
      b2     <- ifelse( as.vector( dhat ) <= eps, tau1, tau1 + omega * as.vector( dhat ) + ( 2 * omega * tau1 ) * t1 )
      b3     <- ifelse( as.vector( dhat ) <= eps, 0, t2 + tau4 * t1 )
      upper  <- alpha2 * b1 + alpha3 * b2 + g * b3
      ksi    <- upper / lower
      tt[[1]]     <- transform( ksi, disobj[[1]], w = as.vector( w ) )
      dhat   <- matrix( tt[[1]]$res, n, m )  ## dhat update
    } else { ## conditionality == "row"
      for ( i in 1:n ) {
        mnc1 <- sum( ps$c1 )
        mnc2 <- sum( ps$c2 )
        nrmw <- sqrt( ps$nrmw2[i] )
        nrmm <- sqrt( ps$nrmm2[i] )
        nrmv <- sqrt( ps$nrmv2[i] )
        sqrtn <- sqrt( n )
        sumc1 <- mnc1 - ps$nrmd2[i]
        sumc2 <- mnc2 - ( ps$nrmv2[i] + omega * ps$nrmm2[i] ) / ps$nrmv2[i]
        g1 = ( sumc1 + ps$nrmd2[i] ) / wpp
        g2 = ( 1 + sumc2 ) * ps$nrmv2[i] + omega * ps$nrmm2[i]
        g3 = n * ps$nrmv2[i]
        g = sqrt ( g1^lambda * g2 / g3 )
        alpha2 <- 0.5 * lambda * sqrt( g2 ) * g1^( lambda/2 - 1 )
        alpha3 <- 0.5 * g1^( lambda / 2 ) / sqrt( g2 )
        fmin <- min( dhat[i, ] )
        tau1 <- ps$wsum[i] / ps$sumw[i]
        tau2   <- 0.5 / eps
        if ( fmin >= eps ) tau2 <- 0.5 / fmin
        tau4   <- tau1 / nrmv
        beta1  <- 1 / wpp
        beta3  <- 1 + sumc2 + omega + 2 * omega * tau1 * tau2
        beta5  <- sqrtn * tau2 * tau4
        lower  <- alpha2 * beta1 + alpha3 * beta3 + g * beta5
        t1     <- tau2 * dhat[i, ] - 0.5
        t2     <- dhat[i, ] / ( 2 * nrmv )
        b1     <- d[i, ] / wpp
        b2     <- ifelse( dhat[i, ] <= eps, tau1 + tau1 * sumc2, tau1 + tau1 * sumc2 + omega * dhat[i, ] + ( 2 * omega * tau1 ) * t1 )
        b3     <- ifelse( dhat[i, ] <= eps, 0, sqrtn * t2 + sqrtn * tau4 * t1 )
        upper  <- alpha2 * b1 + alpha3 * b2 + g * b3
        ksi    <- upper / lower
        tt[[i]] <- transform(ksi, disobj[[i]], w = w[i, ])
        dhat[i, ] <- tt[[i]]$res  ## dhat update
        ps$c1[i] <- ps$nrmd2[i] <- sum( w[i, ] * ( dhat[i, ] - d[i, ] )^2 )
        ps$wsum[i] <- sum( w[i, ] * dhat[i, ] )
        ps$nrmw2[i] <- sum( w[i, ] * dhat[i, ]^2 )
        ps$nrmm2[i] <- ps$wsum[i]^2 / ps$sumw[i]
        ps$nrmv2[i] <- ps$nrmw2[i] - ps$nrmm2[i]
        ps$c2[i] <- ( ps$nrmv2[i] + omega * ps$nrmm2[i] ) / ps$nrmv2[i]
      }
    } # psychometrika target

    ps   <- pstress( dhat, d, w, omega, lambda, wpp )   #pstress value
    lnew <- ps$pstress

    if (verbose) cat("Iteration: ", formatC(itel, digits=6, width=6),
                     "   Stress:",  formatC(lnew, digits=6, width=12, format="f"),
                     "   Dif:",     formatC(lold - lnew, digits=6,width=12, format="f"),
                     "\n")

    #if ( ( (lold-lnew) < eps & itel > 1) || (itel==itmax)) break()
    if ( itel == itmax ) break() 
    if ( 2 * ( lold - lnew ) <= eps * ( lold + lnew + 1e-15 ) ) break()
    
    lold <- lnew                       #update stress
    itel <- itel+1
  }
  #-------------------- end majorization --------------------------

  # point stress
  resmat <- as.matrix(d - dhat)^2    #point stress
  spp.col <- colSums(w * resmat, na.rm = TRUE) / colSums(w, na.rm = TRUE)
  spp.col <- spp.col/sum(spp.col)*100
  spp.row <- rowSums(w * resmat, na.rm = TRUE) / rowSums(w, na.rm = TRUE)
  spp.row <- spp.row/sum(spp.row)*100


  # Add names
  colnames(y) <- colnames(x) <- paste("D",1:(dim(y)[2]),sep="")
  names(spp.row) <- rownames(x) <- rownames(diss) <- rownames(d) <- rnames
  names(spp.col) <- rownames(y) <- colnames(delta)

  # final dilation
  # first re-scale distances (unconditional)
  ssqd <- sum( w * d^2 )
  # if ( ssqd != 0 ) {
    scale <- sqrt( sum( w ) / ssqd )
    d <- scale * d
    x <- scale * x
    y <- scale * y
    dhat <- scale * dhat
  # }
  # then optimally adapt d-hats depending on conditionality
  #if ( conditionality == "matrix" ) {
    # fwd <- sum( dhat * w * d )
    # dwd <- sum( dhat * w * dhat )
    # if (dwd != 0) ( fwd / dwd ) * dhat
  #} 
  #else { ## conditionality == "row"
  #  fwd <- rowSums( dhat * w * d)
  #  dwd <- rowSums( d * w * d)
  #  for ( i in 1:n ) dhat[i] <- ifelse( fwd[i] == 0, dhat[i], ( dwd[i] / fwd[i] ) * dhat[i] )
  #}

  if (itel == itmax) warning("Iteration limit reached! Increase itmax argument!")

  ## stress normalization
  lnew <- sqrt(sum(w*(dhat - d)^2, na.rm = TRUE)/wpp)
  lnew <- sqrt(1 - sum(w * dhat * d, na.rm = TRUE)^2 / 
                  (sum(w * dhat^2, na.rm = TRUE) * sum(w * d^2, na.rm = TRUE)) )

  ## congruence coefficients
  diss0 <- diss
  diss0[is.na(diss0)] <- 0
  congnum <- diag(diss0 %*% t(d))
  congdenom <- sqrt(diag(diss0 %*% t(diss0)) * diag(d %*% t(d)))
  congvec <- congnum/congdenom

  ##
  iord <- vector(mode = "list", length = length(tt))
  if (conditionality == "matrix"){
    iord[[1]] <- tt[[1]]$iord.prim
  } else {
    iord <- vector(mode = "list", n)
    for (i in 1:n) iord[[i]] <- tt[[i]]$iord.prim
  }

  #return configuration distances, row and column configurations, stress
  result <- list(obsdiss = diss, confdist = d, dhat = dhat, iord = iord, conf.row = x, conf.col = y, stress = lnew,
                 pstress = ps$pstress, spp.row = spp.row, spp.col = spp.col, congvec = congvec, weightmat = w,
                 ndim = p, model = "Rectangular smacof", niter = itel, nind = n, nobj = m, trans = trans, conditionality = conditionality, call = match.call())
  class(result) <- c("smacofR")
  result
}
prefscal <- smacofRect