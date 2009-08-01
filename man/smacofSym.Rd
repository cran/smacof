\name{smacofSym}
\alias{smacofSym}

\title{Symmetric smacof}
\description{
Basic smacof on symmetric dissimilarity matrix 
}
\usage{
smacofSym(delta, ndim = 2, weightmat = NULL, init = NULL, metric = TRUE, ties = "primary", verbose = FALSE, relax = 1, modulus = 1, itmax = 1000, eps = 1e-06)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{delta}{Either a symmetric dissimilarity matrix or an object of class \code{"dist"}}
  \item{ndim}{Number of dimensions}
  \item{weightmat}{Optional matrix with dissimilarity weights}
  \item{init}{Matrix with starting values for configurations (optional)}
  \item{metric}{If \code{FALSE} non-metric MDS is performed}
  \item{ties}{Tie specification for non-metric MDS only: \code{"primary"}, \code{"secondary"}, or \code{"tertiary"}}
  \item{verbose}{If \code{TRUE}, intermediate stress is printed out}
  \item{relax}{Relaxed smacof update}
   \item{modulus}{Number of smacof iterations per monotone regression call}
  \item{itmax}{Maximum number of iterations}
  \item{eps}{Convergence criterion}
}

\value{
  \item{obsdiss}{Observed dissimilarities, normalized}
  \item{confdiss}{Configuration dissimilarities}
  \item{conf}{Matrix of final configurations}
  \item{stress.m}{stress value for metric MDS}
  \item{stress.nm}{stress value for non-metric MDS (if computed)}
  \item{ndim}{Number of dimensions}
  \item{model}{Type of smacof model}
  \item{niter}{Number of iterations}
  \item{nobj}{Number of objects}
}
\references{de Leeuw, J. \& Mair, P. (2009). Multidimensional scaling using majorization: 
The R package smacof. Journal of Statistical Software, 31(3), 1-30, \url{http://www.jstatsoft.org/v31/i03/} 
}
\author{Jan de Leeuw and Patrick Mair}

\seealso{\code{\link{smacofConstraint}}, \code{\link{smacofRect}}, \code{\link{smacofIndDiff}}, \code{\link{smacofSphere.primal}}, \code{\link{smacofSphere.dual}}}
\examples{

## simple SMACOF solution for kinship data
data(kinshipdelta)
res <- smacofSym(kinshipdelta)
res
summary(res)

## 3D nonmetric SMACOF solution for trading data
data(trading)
res <- smacofSym(trading, ndim = 3, metric = FALSE, ties = "secondary")
res
}

\keyword{models}
