\name{randomstress}
\alias{randomstress}

\title{Stress Calculation for Random Dissimilarities}
\description{Creates random dissimilarity matrices (n objects), fits an MDS, and returns the stress values 
of each MDS fit. 
}
\usage{
randomstress(n, ndim, nrep = 100, type = c("ratio", "interval", "ordinal", "mspline"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{Number of objects}
  \item{ndim}{Number of dimensions for MDS}
  \item{nrep}{Number of random samples}
  \item{type}{MDS type}
}

\details{The random dissimilarities are drawn from a U(0,1) distribution.
}

\value{   
Returns a vector with stress values.
}

\references{Spence I., Ogilvie, J.C. (1973). A table of expected stress values for random rankings in
nonmetric multidimensional scaling. Multivariate Behavioral Research, 8, 511-517.
}

\examples{
## 8 objects, 2 dimensions, interval MDS (50 replications)
stressvec <- randomstress(n = 8, ndim = 2, nrep = 50, type = "interval")
mean(stressvec)
}

\keyword{models}
