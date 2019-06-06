coefOfVar <- function(x, w){            # Compute coefficient of variation
  av <- sum(x * w)/sum(w)
  va <- sum((x - av)^2 * w)/sum(w)
  return(va/av)
}

# psychometrika stress
pstress <- function( dhat, d, w, omega, lambda, wpp, conditionality ){
  if ( conditionality == "matrix" ){
    c1 <- nrmd2 <- sum( w * ( dhat - d )^2 )
    sumw <- sum( w )
    mnc1 <- nrmd2 / wpp                                    ## raw stress
    wsum <- sum( w * dhat )
    nrmw2 <- sum( w * dhat^2 )
    nrmm2 <- wsum^2 / wpp
    nrmv2 <- nrmw2 - nrmm2                                 ## coefficient of variation (weighted d-hats)
    c2 <- mnc2 <- ( nrmv2 + omega * nrmm2 ) / nrmv2        ## penalty term \mu
    g <- sqrt( mnc1^lambda * mnc2 )                        ## p-stress   
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
