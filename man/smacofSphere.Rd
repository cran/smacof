\name{smacofSphere}
\alias{smacofSphere}


\title{Spherical SMACOF}
\description{Dual and primal approach for spherical SMACOF.
}
\usage{
smacofSphere(delta, ndim = 2, type = c("ratio", "interval", "ordinal","mspline"), 
             algorithm = c("dual", "primal"), weightmat = NULL, 
             init = "torgerson", ties = "primary", verbose = FALSE, penalty = 100, 
             relax = FALSE, modulus = 1, itmax = 1000, eps = 1e-6,
             spline.degree = 2, spline.intKnots = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{delta}{Either a symmetric dissimilarity matrix or an object of class \code{dist}}
  \item{ndim}{Number of dimensions}
  \item{type}{MDS type: \code{"interval"}, \code{"ratio"}, or \code{"ordinal"} (nonmetric MDS)}
  \item{algorithm}{Algorithm type (see details)}
  \item{weightmat}{Optional matrix with dissimilarity weights}
  \item{init}{Either \code{"torgerson"} (classical scaling starting solution), \code{"random"} (random configuration), or a user-defined matrix}
  \item{ties}{Tie specification for non-metric MDS only}
  \item{verbose}{If \code{TRUE}, intermediate stress is printed out}
  \item{penalty}{Penalty parameter for dual algorithm (larger 0), see details}
  \item{relax}{If \code{TRUE}, block relaxation is used for majorization (dual algorithm)}
  \item{modulus}{Number of smacof iterations per monotone regression call}
  \item{itmax}{Maximum number of iterations}
  \item{eps}{Convergence criterion}
  \item{spline.degree}{Degree of the spline for \code{"mspline"} MDS type}
  \item{spline.intKnots}{Number of interior knots of the spline for \code{"mspline"} MDS type}
}

\value{
  \item{delta}{Observed dissimilarities}
  \item{obsdiss}{Observed dissimilarities, normalized}
  \item{obsdiss1}{Dual SMACOF: Observed dissimilarities}
  \item{obsdiss2}{Dual SMACOF: Restriction matrix}
  \item{confdist}{Configuration dissimilarities}
  \item{conf}{Matrix with fitted configurations}
  \item{spp}{Stress per point}
  \item{resmat}{Matrix with squared residuals}
  \item{rss}{Residual sum-of-squares}
  \item{stress}{Stress-1 value}
  \item{init}{Starting configurations}
  \item{ndim}{Number of dimensions}
  \item{dummyvec}{Dummy vector of restriction matrix}
  \item{model}{Type of smacof model}
  \item{niter}{Number of iterations}
  \item{nobj}{Number of objects}
}

\details{
For large scale problems it is suggested to use the dual algorithm. Using the penalty parameter (dual algorithm), the user allow for slight point deviations from the circle (the higher the penalty, the stricter the algorithm is in terms of placing points in the sphere, see examples section below).
}

\references{De Leeuw, J. & Mair, P. (2009). Multidimensional scaling using majorization: 
The R package smacof. Journal of Statistical Software, 31(3), 1-30, \url{https://www.jstatsoft.org/v31/i03/}}
\author{Jan de Leeuw and Patrick Mair}

\seealso{\code{\link{smacofRect}}, \code{\link{smacofIndDiff}}, \code{\link{smacofSym}},\code{\link{smacofConstraint}}}
\examples{

## spherical SMACOF solution for trading data
## dual algorithm
res <- smacofSphere(trading, type = "ordinal")  
res
plot(res)

## lower penalty
res <- smacofSphere(trading, penalty = 20, type = "ordinal")  
res
plot(res)

## primal algorithm, interval
res <- smacofSphere(trading, type = "interval", algorithm = "primal")  
res
}

\keyword{multivariate}
