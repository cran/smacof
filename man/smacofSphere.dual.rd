\name{smacofSphere.dual}
\alias{smacofSphere.dual}
\alias{smacofSphere.primal}


\title{Spherical SMACOF}
\description{Dual and primal approach for spherical SMACOF.
}
\usage{
smacofSphere.dual(delta, penalty = 100, ndim = 2, weightmat = NULL, init = NULL, 
metric = TRUE, ties = "primary", verbose = FALSE, relax = FALSE, modulus = 1, itmax = 1000, eps = 1e-6)

smacofSphere.primal(delta, ndim = 2, weightmat = NULL, init = NULL,
metric = TRUE, ties = "primary", verbose = FALSE, modulus = 1, itmax = 1000, eps = 1e-6)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{delta}{Either a symmetric dissimilarity matrix or an object of class \code{dist}}
  \item{penalty}{Penalty parameter for dual algorithm (larger 0)}
  \item{ndim}{Number of dimensions}
  \item{weightmat}{Optional matrix with dissimilarity weights}
  \item{init}{Matrix with starting values for configurations (optional)}
  \item{metric}{If \code{FALSE} non-metric MDS is performed}
  \item{ties}{Tie specification for non-metric MDS only}
  \item{verbose}{If \code{TRUE}, intermediate stress is printed out}
  \item{relax}{If \code{TRUE}, block relaxation is used for majorization}
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
  \item{conf}{Matrix of final configurations}
  \item{spp}{Stress per point}
  \item{stress.m}{stress value for metric MDS}
  \item{stress.nm}{stress value for non-metric MDS (if computed)}
  \item{ndim}{Number of dimensions}
  \item{dummyvec}{Dummy vector of restriction matrix}
  \item{model}{Type of smacof model}
  \item{niter}{Number of iterations}
  \item{nobj}{Number of objects}
}
\references{de Leeuw, J. \& Mair, P. (2008). Multidimensional scaling using majorization: The R package smacof.}
\author{Jan de Leeuw and Patrick Mair}

\seealso{\code{\link{smacofRect}}, \code{\link{smacofIndDiff}}, \code{\link{smacofSym}},\code{\link{smacofConstraint}}}
\examples{

## spherical SMACOF solution for trading data
data(trading)
res <- smacofSphere.dual(trading)
res
summary(res)


}

\keyword{models}
