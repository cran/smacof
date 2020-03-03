smacofRect <- function(delta, ndim = 2, type = c("ratio", "interval", "ordinal", "mspline"),
                       conditionality = "unconditional", lambda = 0.5, omega = 1, 
                       circle = c("none", "row", "column"), weightmat = NULL, init = NULL, 
                       fixed = c("none", "row", "column"), fixed.coord = NULL,
                       ties = c("primary", "secondary"), verbose = FALSE, relax = TRUE, itmax = 10000,  
                       eps = 1e-6, spline.degree = 2, spline.intKnots = 2, parallelize = FALSE)
  
  # init ... either a list of 2 matrices of dimension n \times p and m \times p with 
  # starting values. if NULL, svd is used.
{
  circle <- match.arg(circle, c("none", "row", "column"), several.ok = FALSE)
  type <- match.arg(type, c("ratio", "interval", "ordinal", "mspline"), several.ok = FALSE)
  ties <- match.arg(ties, c("primary","secondary"), several.ok = FALSE)
  conditionality <- match.arg(conditionality, c("unconditional", "matrix", "row"), several.ok = FALSE)
  if (conditionality == "unconditional") conditionality <- "matrix"         ## for backward compatibility
  fixed <- match.arg(fixed, c("none", "row", "column"), several.ok = FALSE)
  diss <- delta
  rnames <- rownames(delta)
  if (is.data.frame(diss)) diss <- as.matrix(diss)
  checkdiss(diss)

  n <- dim(diss)[1]                       #number of individuals
  m <- dim(diss)[2]                       #number of objects
  p <- ndim
  
  ## --- sanity check external unfolding
  if (is.null(fixed.coord) & fixed != "none") stop("fixed.coord is not specified.")
  if (fixed == "row"){
    if (nrow(fixed.coord) != n | ncol(fixed.coord) < p) 
      stop("Matrix of fixed row coordinates specified by fixed.coord has an incorrect size")
  } else if (fixed == "col"){
    if (nrow(fixed.coord) != m | ncol(fixed.coord) < p) 
      stop("Matrix of fixed column coordinates specified by fixed.coord has an incorrect size")
  }  

  ## --- weight init
  if (is.null(weightmat)) {
    w <- matrix(1,n,m)                      #initialize weights (as 1)
    if (any(is.na(diss))) w[is.na(diss)] <- 0  # blank out missings
  } else w <- weightmat
  wpp <- sum(w)

  ## --- Prepare for parallization
  if (parallelize) {
    noCores <- detectCores()/2
    cl <- parallel::makeCluster(noCores)
    #cl <- parallel::makeCluster(4)
    doParallel::registerDoParallel(cl)
  } 
  
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
  ## --- end optimal scaling prep
  
  ## --- conditionality prep
  disobj <- list()
  if (conditionality == "matrix"){
    disobj[[1]] <- transPrep(as.vector(delta), trans = trans, spline.intKnots = spline.intKnots, spline.degree = spline.degree)
    if (trans == "mspline") disobj[[1]]$base <- cbind(rep(1, nrow(disobj[[1]]$base)), disobj[[1]]$base)
    tt <- transform(as.vector(delta), disobj[[1]], w = as.vector(w)) #, normq = wpp
    dhat <- matrix(tt$res, n, m)  ## dhat update
    tt <- vector(mode = "list", length = 1)
  } else {                             ## conditionality == "row"
    dhat <- matrix(0, n, m)
    for (i in 1:n) {
      disobj[[i]] <- transPrep(delta[i, ], trans = trans, spline.intKnots = spline.intKnots, spline.degree = spline.degree)
      if (trans == "mspline") disobj[[i]]$base <- cbind(rep(1, nrow(disobj[[i]]$base)), disobj[[i]]$base)
      dhat[i, ] <- transform(delta[i, ], disobj[[i]], w = w[i, ], normq = sum(w))$res  ## dhat update
    }
    tt <- vector(mode = "list", length = n)
  }
 ## --- end conditionality prep


  itel <- 1
  e <- dhat^2
  e <- -0.5*(e - outer(rowSums(e)/m, colSums(e)/n, "+") + (sum(e)/(n*m)))
  
  ## --------- external unfolding ------------
  if (is.list(init) & fixed == "none") {          ## Initial coordinates given
    x <- init[[1]][, 1:p]                         
    y <- init[[2]][, 1:p]
  } else if (is.list(init) & fixed == "row") {    ## Initial coordinates given + fixed rows
    fixed.coord <- fixed.coord[, 1:p]
    x <- fixed.coord            
    y <- init[[2]][, 1:p]
  } else if (is.list(init) & fixed == "column") {    ## Initial coordinates given + fixed columns
    fixed.coord <- fixed.coord[, 1:p]
    x <- init[[1]][, 1:p]                     
    y <- fixed.coord
  } else if (!is.list(init) & fixed == "none") {  ## No initial coordinates given
    z <- svd(e, nu = p, nv = 0)           #SVD for e (pos. distances)
    x <- z$u                              #starting value for x
    y <- crossprod(e, x)                  #starting value for y
  } else if (!is.list(init) & fixed == "row") {  ## No initial coordinates given + fixed rows
    fixed.coord <- fixed.coord[, 1:p]
    x <- fixed.coord                       #starting value for x
    y <- crossprod(e, x) %*% solve(crossprod(x)) #starting value for y
  } else if (!is.list(init) & fixed == "column") {  ## No initial coordinates given + fixed cols
    fixed.coord <- fixed.coord[, 1:p]
    y <- fixed.coord                       #starting value for y
    x <- solve(crossprod(y), crossprod(y, t(e))) #starting value for x
    x <- t(x)
  }
  ## Rescale x and y by a constant to fit e
  dil <- sum(diag(t(x) %*% (e %*% y))) / sum(diag(crossprod(x) %*% crossprod(y)))
  x   <- (dil^.5/ sum(x^2)^.5 ) * x
  y   <- (dil^.5/ sum(y^2)^.5 ) * y

  
  ## --------- circular restriction 
  if (circle != "none"){
    r <- projCircle(x,y,x,y,circle=circle)
    x <- r$x
    y <- r$y
    wr <- rowSums(w)
    wc <- colSums(w)
    lambda.circ <- 2*max(c(wr,wc))
  }

  d <- distRect(x,y,0)                  #n times m of reproduced diss
  ps   <- pstress(dhat, d, w, omega, lambda, wpp, conditionality)   #pstress value
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
      yunc <- v %*% (yraw + crossprod(ww, xraw/wr))
      xunc <- (xraw + (ww %*% yunc))/wr
      ## x update
      if (fixed != "column") {
        y <- yunc                
      } else {
        y <- (sum(diag(crossprod(fixed.coord, yunc))/sum(fixed.coord^2))) * fixed.coord
      }
      ## y update
      if (fixed != "row") {
        x <- xunc        
      } else {
        x <- (sum(diag(crossprod(fixed.coord, xunc))/sum(fixed.coord^2))) * fixed.coord
      }
      if (relax & itel > 100){
        x <- 2 * x - xold
        y <- 2 * y - yold
      }
      
    } else {
      b  <- w * (1 - (!d.is.0) * dhat / (d + d.is.0))  #B matrix
      br <- rowSums(b)                   #rows B
      bc <- colSums(b)                   #columns W
      xunc <- x - outer(br, rep(1/lambda.circ, p), "*") * x + b %*% (y/lambda.circ)
      yunc <- y - outer(bc, rep(1/lambda.circ, p), "*") * y + t(b) %*% (x/lambda.circ)
      r <- projCircle(xunc, yunc, x, y, circle = circle)
      if (fixed != "row") {
        x <- r$x
      } else {
        x <- (sum(diag(crossprod(fixed.coord, xunc))/sum(fixed.coord^2))) * fixed.coord
      }
      if (fixed != "column") {
        y <- r$y
      } else {
        y <- (sum(diag(crossprod(fixed.coord, yunc))/sum(fixed.coord^2))) * fixed.coord
      }
      
    }

    d <- distRect(x, y, 0)               #compute distances (update)
    dhat.old <- dhat

    
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
      #b2     <- ifelse( as.vector( dhat ) <= eps, tau1, tau1 + omega * as.vector( dhat ) + ( 2 * omega * tau1 ) * t1 )
      b2     <- tau1 +  (as.vector(dhat) > eps) * omega * (dhat + ( 2 * tau1 ) * t1)
      #b3     <- ifelse( as.vector( dhat ) <= eps, 0, t2 + tau4 * t1 )
      b3     <- (as.vector(dhat) > eps) * (t2 +  tau4 * t1)
      upper  <- alpha2 * b1 + alpha3 * b2 + g * b3
      ksi    <- upper / lower
      ##tt[[1]]     <- transform( ksi, disobj[[1]], w = as.vector( w ) )
      tt[[1]]     <- smacof::transform( ksi, disobj[[1]], w = as.vector( w ) )
      dhat   <- matrix( tt[[1]]$res, n, m )  ## dhat update
    } else {                      ## conditionality == "row"
      if (parallelize){
        ii <- c(0, floor(n * (1:noCores)/noCores))
        s <- 2:(noCores + 1)
        tt <- foreach(s = s, .combine = c, .verbose = FALSE, 
                      .maxcombine = noCores, .multicombine = TRUE, export = ps) %dopar% {
          dhatRowUpdateBlock(ps, s, ii, n, omega, lambda, wpp, eps, dhat, d, w, disobj)
        }
        for (i in 1:n) {
          dhat[i, ] <- tt[[i]]$res
        }
      } else {
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
          #b2     <- ifelse( dhat[i, ] <= eps, tau1 + tau1 * sumc2, tau1 + tau1 * sumc2 + omega * dhat[i, ] + ( 2 * omega * tau1 ) * t1 )
          b2     <- (tau1 + tau1 * sumc2) + (dhat[i, ] > eps) * omega * (dhat[i, ] + ( 2 * tau1 ) * t1)
          #b2     <- ifelse( dhat[i, ] <= eps, tau1 + tau1 * sumc2, tau1 + tau1 * sumc2 + omega * dhat[i, ] + ( 2 * omega * tau1 ) * t1 )
          #b3     <- ifelse( dhat[i, ] <= eps, 0, sqrtn * (t2 +  tau4 * t1) )
          b3     <- (dhat[i, ] > eps) * sqrtn * (t2 +  tau4 * t1)
          upper  <- alpha2 * b1 + alpha3 * b2 + g * b3
          ksi    <- upper / lower
          tt[[i]] <- transform(ksi, disobj[[i]], w = w[i, ])
          dhat[i, ] <- tt[[i]]$res  ## dhat update
          # ps$c1[i] <- ps$nrmd2[i] <- sum( w[i, ] * ( dhat[i, ] - d[i, ] )^2 )
          # ps$wsum[i] <- sum( w[i, ] * dhat[i, ] )
          # ps$nrmw2[i] <- sum( w[i, ] * dhat[i, ]^2 )
          # ps$nrmm2[i] <- ps$wsum[i]^2 / ps$sumw[i]
          # ps$nrmv2[i] <- ps$nrmw2[i] - ps$nrmm2[i]
          # ps$c2[i] <- ( ps$nrmv2[i] + omega * ps$nrmm2[i] ) / ps$nrmv2[i]
        }
      }
    }
    
    ps   <- pstress( dhat, d, w, omega, lambda, wpp, conditionality )   #pstress value
    lnew <- ps$pstress

    if (verbose) cat("Iteration: ", formatC(itel, digits=6, width=6),
                     "  P-Stress:",  formatC(lnew, digits=6, width=12, format="f"),
                     "   Dif:",     formatC(lold - lnew, digits=6,width=12, format="f"),
                     "\n")

    #if ( ( (lold-lnew) < eps & itel > 1) || (itel==itmax)) break()
    if ( itel == itmax ) break() 
    if ( 2 * ( lold - lnew ) <= eps * ( lold + lnew + 1e-15 ) ) break()
    
    lold <- lnew                       #update stress
    itel <- itel+1
  }
  #-------------------- end majorization --------------------------

  ## --- Finish for parallization
  if (parallelize) parallel::stopCluster(cl)
  
  
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
  scale <- sqrt( sum( w ) / ssqd )
  d <- scale * d
  x <- scale * x
  y <- scale * y
  dhat <- scale * dhat

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
  return(result)
}
