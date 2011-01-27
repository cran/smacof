\name{smacofRect}
\alias{smacofRect}

\title{Rectangular smacof}
\description{Variant of smacof for rectangular matrices (typically ratings, preferences) which is also known as metric unfolding.
}
\usage{
smacofRect(delta, ndim = 2, weightmat = NULL, init = NULL, verbose = FALSE,
           itmax = 1000, reg = 1e-6, eps = 1e-6)
}

\arguments{
  \item{delta}{Data frame or matrix of preferences, ratings, dissimilarities.}
  \item{ndim}{Number of dimensions}
  \item{weightmat}{Optional matrix with dissimilarity weights}
  \item{init}{Matrix with starting values for configurations (optional)}
  \item{verbose}{If \code{TRUE}, intermediate stress is printed out}
  \item{itmax}{Maximum number of iterations}
  \item{reg}{Regularization factor, prevents distances from beeing 0}
  \item{eps}{Convergence criterion}
}
\details{Creates an object of class \code{smacofR}. 

}

\value{
  \item{obsdiss}{Observed dissimilarities, corresponds to delta}
  \item{confdiss}{Configuration dissimilarities}
  \item{conf.row}{Matrix of final row configurations}
  \item{conf.col}{Matrix of final column configurations}
  \item{stress}{Final stress value}
  \item{spp.row}{Stress per point, rows}
  \item{spp.col}{Stress per point, columns}
  \item{ndim}{Number of dimensions}
  \item{model}{Type of smacof model}
  \item{niter}{Number of iterations}
  \item{nind}{Number of individuals (rows)}
  \item{nobj}{Number of objects (columns)}
}
\references{de Leeuw, J. \& Mair, P. (2009). Multidimensional scaling using majorization: 
The R package smacof. Journal of Statistical Software, 31(3), 1-30, \url{http://www.jstatsoft.org/v31/i03/} 
}
\author{Jan de Leeuw and Patrick Mair}

\seealso{\code{\link{smacofConstraint}}, \code{\link{smacofSym}}, \code{\link{smacofIndDiff}}, \code{\link{smacofSphere.primal}}, \code{\link{smacofSphere.dual}}}

\examples{
data(breakfast)
res <- smacofRect(breakfast)
res
summary(res)
}

\keyword{models}
