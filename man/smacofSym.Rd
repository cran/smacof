\name{smacofSym}
\alias{smacofSym}
\alias{mds}

\title{Symmetric smacof}
\description{
Multidimensional scaling on a symmetric dissimilarity matrix using SMACOF. 
}
\usage{
smacofSym(delta, ndim = 2, type = c("ratio", "interval", "ordinal", "mspline"), 
          weightmat = NULL, init = "torgerson", ties = "primary", principal = FALSE, 
          verbose = FALSE, relax = FALSE, modulus = 1, itmax = 1000, eps = 1e-06, 
          spline.degree = 2, spline.intKnots = 2)

mds(delta, ndim = 2, type = c("ratio", "interval", "ordinal", "mspline"), 
    weightmat = NULL, init = "torgerson", ties = "primary", principal = FALSE, 
    verbose = FALSE, relax = FALSE, modulus = 1, itmax = 1000, eps = 1e-06, 
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
  \item{principal}{If \code{TRUE}, principal axis transformation is applied to the final configuration}
  \item{verbose}{If \code{TRUE}, intermediate stress is printed out}
  \item{relax}{If \code{TRUE}, block relaxation is used for majorization}
  \item{modulus}{Number of smacof iterations per monotone regression call}
  \item{itmax}{Maximum number of iterations}
  \item{eps}{Convergence criterion}
  \item{spline.degree}{Degree of the spline for \code{"mspline"} MDS type}
  \item{spline.intKnots}{Number of interior knots of the spline for \code{"mspline"} MDS type}
}
\details{
The function \code{mds()} is a wrapper function and can be used instead of \code{smacofSym()}. It reports the Stress-1 value (normalized). The main output are the coordinates in the low-dimensional space (configuration; \code{conf}; see also \code{\link{plot.smacof}}). 

Four types of MDS can be fitted: ratio MDS (no dissimilarity transformation), interval MDS (linear transformation), ordinal MDS (ordinal transformation with various options for handling ties), and spline MDS (monotone spline transformation). Shepard plots in \code{\link{plot.smacof}} give insight into this transformation.

Setting \code{principal = TRUE} is useful for interpretatbility of the dimensions, or to check hypotheses about the dimensions.

In case of missing input dissimilarities, the \code{weightmat} is computed internally so that missings are blanked out during optimization. 
}

\value{
  \item{delta}{Observed dissimilarities, not normalized}
  \item{dhat}{Disparities (transformed proximities, approximated distances, d-hats)}
  \item{confdist}{Configuration distances}
  \item{conf}{Matrix of fitted configurations}
  \item{stress}{Stress-1 value}
  \item{spp}{Stress per point (stress contribution of each point on a percentage scale)}
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
\references{
De Leeuw, J. & Mair, P. (2009). Multidimensional scaling using majorization: 
The R package smacof. Journal of Statistical Software, 31(3), 1-30, \doi{10.18637/jss.v031.i03} 

Mair, P, Groenen, P. J. F., De Leeuw, J. (2022). More on multidimensional scaling in R: smacof version 2. Journal of Statistical Software, 102(10), 1-47. \doi{10.18637/jss.v102.i10}

Borg, I., & Groenen, P. J. F. (2005). Modern Multidimensional Scaling (2nd ed.). Springer.

Borg, I., Groenen, P. J. F., & Mair, P. (2018). Applied Multidimensional Scaling and Unfolding (2nd ed.). Springer. 
}

\author{Jan de Leeuw, Patrick Mair, and Patrick Groenen}

\seealso{\code{\link{smacofConstraint}}, \code{\link{smacofRect}}, \code{\link{smacofIndDiff}}, \code{\link{smacofSphere}}, \code{\link{plot.smacof}}}
\examples{

## simple SMACOF solution (interval MDS) for kinship data
res <- mds(kinshipdelta, type = "interval")
res
summary(res)
plot(res)
plot(res, type = "p", label.conf = list(label = TRUE, col = "darkgray"), pch = 25, col = "red")

## ratio MDS, random starts
set.seed(123)
res <- mds(kinshipdelta, init = "random")
res

## 3D ordinal SMACOF solution for trading data (secondary approach to ties)
data(trading)
res <- mds(trading, ndim = 3, type = "ordinal", ties = "secondary")
res

## spline MDS 
delta <- sim2diss(cor(PVQ40agg))
res <- mds(delta, type = "mspline", spline.degree = 3, spline.intKnots = 4)
res
plot(res, "Shepard")
}

\keyword{multivariate}
