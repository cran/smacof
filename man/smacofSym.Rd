\name{smacofSym}
\alias{smacofSym}

\title{Symmetric smacof}
\description{
Multidimensional scaling (stress minimization: SMACOF) on symmetric dissimilarity matrix. 
}
\usage{
smacofSym(delta, ndim = 2, weightmat = NULL, init = NULL, metric = TRUE, ties = "primary", verbose = FALSE, relax = FALSE, modulus = 1, itmax = 1000, eps = 1e-06)

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
  \item{relax}{If \code{TRUE}, block relaxation is used for majorization}
  \item{modulus}{Number of smacof iterations per monotone regression call}
  \item{itmax}{Maximum number of iterations}
  \item{eps}{Convergence criterion}
}
\details{This is the simplest MDS-SMACOF version of the package. It solves the stress target function for symmetric dissimiliby means of the majorization approach (SMACOF). The main output are the coordinates in the low-dimensional space (configurations; \code{conf}).
The function also returns the point stress, i.e. the larger the contribution of a point to the total stress, the worse the fit (see also \code{\link{plot.smacof}}. 
}
\value{
  \item{delta}{Observed dissimilarities}
  \item{obsdiss}{Observed dissimilarities, normalized}
  \item{confdiss}{Configuration dissimilarities}
  \item{conf}{Matrix of final configurations}
  \item{stress.m}{Stress value for metric MDS}
  \item{stress.nm}{Stress value for non-metric MDS (if computed)}
  \item{spp}{Stress per point}
  \item{ndim}{Number of dimensions}
  \item{model}{Type of smacof model}
  \item{niter}{Number of iterations}
  \item{nobj}{Number of objects}
}
\references{de Leeuw, J. & Mair, P. (2009). Multidimensional scaling using majorization: 
The R package smacof. Journal of Statistical Software, 31(3), 1-30, \url{http://www.jstatsoft.org/v31/i03/} 
}
\author{Jan de Leeuw and Patrick Mair}

\seealso{\code{\link{smacofConstraint}}, \code{\link{smacofRect}}, \code{\link{smacofIndDiff}}, \code{\link{smacofSphere.primal}}, \code{\link{smacofSphere.dual}}, \code{\link{plot.smacof}}}
\examples{

## simple SMACOF solution for kinship data
data(kinshipdelta)
res <- smacofSym(kinshipdelta)
res
summary(res)
plot(res)
plot(res, type = "p", label.conf = list(TRUE, 3, "darkgray"), pch = 25, col = "red")

## 3D nonmetric SMACOF solution for trading data
data(trading)
res <- smacofSym(trading, ndim = 3, metric = FALSE, ties = "secondary")
res
}

\keyword{models}
