`initWeights` <-
function(diss) {
  if (!is.list(diss)) {
	n<-attr(diss,"Size")
	return(as.dist(matrix(1,n,n)))
  }
  n <- attr(diss[[1]],"Size"); m<-length(diss)
  return(repList(as.dist(matrix(1,n,n)),m))
}

