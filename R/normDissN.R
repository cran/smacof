`normDissN` <-
function(diss,wghts,m)
{
  N <- length(diss)*m
  dissnorm <- diss/sqrt(sum(wghts*diss^2))*sqrt(N)
  return(dissnorm)
}

