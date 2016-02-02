\name{smacofSym}
\alias{smacofSym}
\alias{mds}

\title{Symmetric smacof}
\description{
Multidimensional scaling on a symmetric dissimilarity matrix using SMACOF. 
}
\usage{
smacofSym(delta, ndim = 2, type = c("ratio", "interval", "ordinal", "mspline"), 
          weightmat = NULL, init = "torgerson", ties = "primary", verbose = FALSE, 
          relax = FALSE, modulus = 1, itmax = 1000, eps = 1e-06, 
          spline.degree = 2, spline.intKnots = 2)

mds(delta, ndim = 2, type = c("ratio", "interval", "ordinal", "mspline"), 
    weightmat = NULL, init = "torgerson", ties = "primary", verbose = FALSE, 
    relax = FALSE, modulus = 1, itmax = 1000, eps = 1e-06, 
    spline.degree = 2, spline.intKnots = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{delta}{Either a symmetric dissimilarity matrix or an object of class \code{"dist"}}
  \item{ndim}{Number of dimensions}
  \item{weightmat}{Optional matrix with dissimilarity weights}
  \item{init}{Either \code{"torgerson"} (classical scaling starting solution), \code{"random"} (random configuration), or a user-defined matrix}
  \item{type}{MDS type: \code{"interval"}, \code{"ratio"}, \code{"ordinal"} (nonmetric MDS), or \code{"mspline"}}
  \item{ties}{Tie specification (ordinal MDS only): \code{"primary"}, \code{"secondary"}, or \code{"tertiary"}}
  \item{verbose}{If \code{TRUE}, intermediate stress is printed out}
  \item{relax}{If \code{TRUE}, block relaxation is used for majorization}
  \item{modulus}{Number of smacof iterations per monotone regression call}
  \item{itmax}{Maximum number of iterations}
  \item{eps}{Convergence criterion}
  \item{spline.degree}{Degree of the spline for \code{"mspline"} MDS type}
  \item{spline.intKnots}{Number of interior knots of the spline for \code{"mspline"} MDS type}
}
\details{This is the simplest MDS-SMACOF version of the package. It solves the stress target function for symmetric dissimiliarities by means of the majorization approach (SMACOF) and reports the Stress-1 value (normalized). The main output are the coordinates in the low-dimensional space (configurations; \code{conf}). 

The function \code{mds()} is a wrapper function and can be used instead of \code{smacofSym()}

This function allows for fitting three basic types of MDS: ratio MDS (default), interval MDS (polynomial transformation), and ordinal MDS (aka nonmetric MDS). It also returns the point stress, i.e. the larger the contribution of a point to the total stress, the worse the fit (see also \code{\link{plot.smacof}}. 
}
\value{
  \item{delta}{Observed dissimilarities, not normalized}
  \item{dhat}{Disparities (transformed proximities, approximated distances, d-hats)}
  \item{confdiss}{Configuration distances}
  \item{conf}{Matrix of fitted configurations}
  \item{stress}{Stress-1 value}
  \item{spp}{Stress per point (stress contribution in percentages)}
  \item{resmat}{Matrix with squared residuals}
  \item{rss}{Residual sum-of-squares}
  \item{weightmat}{Weight matrix}
  \item{ndim}{Number of dimensions}
  \item{init}{Starting configuration}
  \item{model}{Name of smacof model}
  \item{niter}{Number of iterations}
  \item{nobj}{Number of objects}
  \item{type}{Type of MDS model}
}
\references{De Leeuw, J. & Mair, P. (2009). Multidimensional scaling using majorization: 
The R package smacof. Journal of Statistical Software, 31(3), 1-30, \url{http://www.jstatsoft.org/v31/i03/} 

Borg, I., & Groenen, P. J. F. (2005). Modern Multidimensional Scaling (2nd ed.). Springer.

Borg, I., Groenen, P. J. F., & Mair, P. (2013). Applied Multidimensional Scaling. Springer. 
}
\author{Jan de Leeuw and Patrick Mair}

\seealso{\code{\link{smacofConstraint}}, \code{\link{smacofRect}}, \code{\link{smacofIndDiff}}, \code{\link{smacofSphere}}, \code{\link{plot.smacof}}}
\examples{

## simple SMACOF solution for kinship data
res <- mds(kinshipdelta)
res
summary(res)
plot(res)
plot(res, type = "p", label.conf = list(label = TRUE, col = "darkgray"), pch = 25, col = "red")

## interval MDS, random starts
res <- smacofSym(kinshipdelta, type = "interval", init = "random")
res

## 3D ordinal SMACOF solution for trading data (secondary approach to ties)
data(trading)
res <- smacofSym(trading, ndim = 3, type = "ordinal", ties = "secondary")
res
}

\keyword{models}
