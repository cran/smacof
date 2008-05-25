`torgerson` <-
function(diss, p=2)
{
#diss ... dissimilarity matrix
#p ... number of dimensions

#------------------- begin subroutines -------------------
#Torgerson's double centering
  doubleCenter <- function(x) {
    n <- dim(x)[1]
    m <- dim(x)[2]
    s <- sum(x)/(n*m)
    xr <- rowSums(x)/m
    xc <- colSums(x)/n
    return((x-outer(xr,xc,"+"))+s)
  }
#-------------------- end subroutines --------------------
  z <- eigen(-doubleCenter(as.matrix(diss)^2)/2,symmetric=TRUE)
  v <- pmax(z$values,0)
  return(z$vectors[,1:p]%*%diag(sqrt(v[1:p])))
}

