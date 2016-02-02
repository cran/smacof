\name{smacofSphere}
\alias{smacofSphere}


\title{Spherical SMACOF}
\description{Dual and primal approach for spherical SMACOF.
}
\usage{
smacofSphere(delta, algorithm = c("dual", "primal"), ndim = 2, 
             type = c("ratio", "interval", "ordinal"), weightmat = NULL, 
             init = "torgerson", ties = "primary", verbose = FALSE, penalty = 100, 
             relax = FALSE, modulus = 1, itmax = 1000, eps = 1e-6)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{delta}{Either a symmetric dissimilarity matrix or an object of class \code{dist}}
  \item{algorithm}{Algorithm type}
  \item{penalty}{Penalty parameter for dual algorithm (larger 0)}
  \item{ndim}{Number of dimensions}
  \item{weightmat}{Optional matrix with dissimilarity weights}
  \item{init}{Either \code{"torgerson"} (classical scaling starting solution), \code{"random"} (random configuration), or a user-defined matrix}
  \item{type}{MDS type: \code{"interval"}, \code{"ratio"}, or \code{"ordinal"} (nonmetric MDS)}
  \item{ties}{Tie specification for non-metric MDS only}
  \item{verbose}{If \code{TRUE}, intermediate stress is printed out}
  \item{relax}{If \code{TRUE}, block relaxation is used for majorization (dual algorith)}
  \item{modulus}{Number of smacof iterations per monotone regression call}
  \item{itmax}{Maximum number of iterations}
  \item{eps}{Convergence criterion}
}

\value{
  \item{delta}{Observed dissimilarities}
  \item{obsdiss}{Observed dissimilarities, normalized}
  \item{obsdiss1}{Dual SMACOF: Observed dissimilarities}
  \item{obsdiss2}{Dual SMACOF: Restriction matrix}
  \item{confdiss}{Configuration dissimilarities}
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
For large scale problems it is suggested to use the dual algorithm. The primal algorithm gets pretty slow, especially for 
ordinal scaling.
}

\references{De Leeuw, J. & Mair, P. (2009). Multidimensional scaling using majorization: 
The R package smacof. Journal of Statistical Software, 31(3), 1-30, \url{http://www.jstatsoft.org/v31/i03/}}
\author{Jan de Leeuw and Patrick Mair}

\seealso{\code{\link{smacofRect}}, \code{\link{smacofIndDiff}}, \code{\link{smacofSym}},\code{\link{smacofConstraint}}}
\examples{

## spherical SMACOF solution for trading data
data(trading)
res <- smacofSphere(trading)
res
summary(res)


}

\keyword{models}
