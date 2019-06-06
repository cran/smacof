dhatRowUpdateBlock <- function(ps, s, ii, n, omega, lambda, wpp, eps, dhat, d, w, disobj){
  tt <- vector(mode = "list", length = ii[s] - ii[s - 1])
  for (i in (ii[s-1] + 1):ii[s]){
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
    #b3     <- ifelse( dhat[i, ] <= eps, 0, sqrtn * t2 + sqrtn * tau4 * t1 )
    b3     <- (dhat[i, ] > eps) * sqrtn * (t2 +  tau4 * t1)
    upper  <- alpha2 * b1 + alpha3 * b2 + g * b3
    ksi    <- upper / lower
    tt[[i - ii[s-1]]]  <- transform(ksi, disobj[[i]], w = w[i, ])
 #   <- dhatRowUpdate(ps, i, n, omega, lambda, wpp, eps, dhat = dhat[i, ], d_i = d[i, ], w_i = w[i, ], disobj_i = disobj[[i]])
  }
  return(tt)
}
# dhatRowUpdate <- function(ps, i, n, omega, lambda, wpp, eps, dhat_i, d_i, w_i, disobj_i){
#   mnc1 <- sum( ps$c1 )
#   mnc2 <- sum( ps$c2 )
#   nrmw <- sqrt( ps$nrmw2[i] )
#   nrmm <- sqrt( ps$nrmm2[i] )
#   nrmv <- sqrt( ps$nrmv2[i] )
#   sqrtn <- sqrt( n )
#   sumc1 <- mnc1 - ps$nrmd2[i]
#   sumc2 <- mnc2 - ( ps$nrmv2[i] + omega * ps$nrmm2[i] ) / ps$nrmv2[i]
#   g1 = ( sumc1 + ps$nrmd2[i] ) / wpp
#   g2 = ( 1 + sumc2 ) * ps$nrmv2[i] + omega * ps$nrmm2[i]
#   g3 = n * ps$nrmv2[i]
#   g = sqrt ( g1^lambda * g2 / g3 )
#   alpha2 <- 0.5 * lambda * sqrt( g2 ) * g1^( lambda/2 - 1 )
#   alpha3 <- 0.5 * g1^( lambda / 2 ) / sqrt( g2 )
#   fmin <- min( dhat_i )
#   tau1 <- ps$wsum[i] / ps$sumw[i]
#   tau2   <- 0.5 / eps
#   if ( fmin >= eps ) tau2 <- 0.5 / fmin
#   tau4   <- tau1 / nrmv
#   beta1  <- 1 / wpp
#   beta3  <- 1 + sumc2 + omega + 2 * omega * tau1 * tau2
#   beta5  <- sqrtn * tau2 * tau4
#   lower  <- alpha2 * beta1 + alpha3 * beta3 + g * beta5
#   t1     <- tau2 * dhat_i - 0.5
#   t2     <- dhat_i / ( 2 * nrmv )
#   b1     <- d_i / wpp
#   b2     <- ifelse( dhat_i <= eps, tau1 + tau1 * sumc2, tau1 + tau1 * sumc2 + omega * dhat_i + ( 2 * omega * tau1 ) * t1 )
#   b3     <- ifelse( dhat_i <= eps, 0, sqrtn * t2 + sqrtn * tau4 * t1 )
#   upper  <- alpha2 * b1 + alpha3 * b2 + g * b3
#   ksi    <- upper / lower
#   tt <- transform(ksi, disobj_i, w = w_i)
#   return(tt)
#   #dhat[i, ] <- tt[[i]]$res  ## dhat update
#   #list(dhat_i, tt_i)
# }