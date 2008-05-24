`bmat` <-
function(diss, wgths, d)
{
  b <- as.matrix((wgths*diss)/(d+ifelse(d==0,1,0)))
  r <- rowSums(b) 
  return(diag(r)-b)
}

