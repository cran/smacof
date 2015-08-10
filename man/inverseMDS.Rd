\name{inverseMDS}
\alias{inverseMDS}

\title{Inverse MDS}
\description{
Performs inverse MDS. That is, given a configuration matrix, it computes a dissimilarity matrix.
}
\usage{
inverseMDS(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Configuration matrix}
}

\value{
Returns a list with dissimilarity matrices. 
}

\details{
Note that this is a very basic implementation. It can be that some dissimilarities are negative. Also, the weights 
are fixed to a value of 1. 
}

\references{
De Leeuw J., Groenen P.J.F. (1997). Inverse multidimensional scaling. Journal of Classification, 12, 3-21.

De Leeuw, J. (2012). On inverse multidimensional scaling. Statistics preprint series, UCLA. \url{gifi.stat.ucla.edu/janspubs/2012/
notes/deleeuw_U_12b.pdf}.
}

\seealso{\code{\link{mds}}}

\examples{
## MDS and inverse MDS on kinship data subset
D <- as.matrix(kinshipdelta)[1:6, 1:6]
fit <- mds(D)                       ## MDS D --> conf
ifit <- inverseMDS(fit$conf)        ## inverse MDS conf --> D

## MDS fit on dissimilarity matrix obtained by inverse MDS
op <- par(mfrow = c(3,3))           
plot(fit, main = "Original MDS")
for (i in 1:length(ifit)) {
  fit.i <- mds(ifit[[i]])             
  plot(fit.i, main = paste0("Inverse MDS (",i, ")"))
}
par(op)

}

\keyword{models}
