\name{smacofIndDiff}
\alias{smacofIndDiff}

\title{smacof for Individual Differences}
\description{Performs smacof for individual differences also known as Three-Way smacof on a list of dissimilarity matrices. Various restrictions decompositions and restrictions on the weight matrix are provided.
}

\usage{
smacofIndDiff(delta, ndim = 2, weightmat = NULL, init = NULL, metric = TRUE,
              ties = "primary", constraint = NULL, verbose = FALSE, modulus = 1,
              itmax = 1000, eps = 1e-6)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{delta}{A list of dissimilarity matrices or a list objects of class \code{dist}}
  \item{ndim}{Number of dimensions}
  \item{weightmat}{Optional matrix with dissimilarity weights}
  \item{init}{Matrix with starting values for configurations (optional)}
  \item{metric}{If \code{FALSE} non-metric MDS is performed}
  \item{ties}{Tie specification for non-metric MDS}
  \item{constraint}{Either \code{NULL}, \code{"idioscal"}, \code{"diagonal"}, or \code{"identity"} (see details)}
  \item{verbose}{If \code{TRUE}, intermediate stress is printed out}
  \item{modulus}{Number of smacof iterations per monotone regression call}
  \item{itmax}{Maximum number of iterations}
  \item{eps}{Convergence criterion}
}
\details{If the constraint is \code{NULL}, INDSCAL is performed with identity configuration weight matrices.
An additional restriction can be imposed with \code{"identity"} which restricts the configurations across individuals/replications/ways to be equal.
More unrestricted models are \code{"diagonal"} which restricts only the configuration weight matrices to be diagonal and \code{"idioscal"} is unrestricted.
}

\value{
  \item{delta}{Observed dissimilarities}
  \item{obsdiss}{List of observed dissimilarities, normalized}
  \item{confdiss}{List of configuration dissimilarities}
  \item{conf}{List of matrices of final configurations}
  \item{gspace}{Joint configurations aka group stimulus space}
  \item{cweights}{Configuration weights}
  \item{stress.m}{stress value for metric MDS}
  \item{stress.nm}{stress value for non-metric MDS (if computed)}
  \item{stress.co}{Constrained stress value}
  \item{spp}{Stress per point}
  \item{ndim}{Number of dimensions}
  \item{model}{Type of smacof model}
  \item{niter}{Number of iterations}
  \item{nobj}{Number of objects}

}
\references{de Leeuw, J. \& Mair, P. (2009). Multidimensional scaling using majorization: 
The R package smacof. Journal of Statistical Software, 31(3), 1-30, \url{http://www.jstatsoft.org/v31/i03/} 
}
\author{Jan de Leeuw and Patrick Mair}

\seealso{\code{\link{smacofConstraint}}, \code{\link{smacofSym}}, \code{\link{smacofRect}}, \code{\link{smacofSphere.primal}}, \code{\link{smacofSphere.dual}}}

\examples{
data(perception)
res <- smacofIndDiff(perception)
res
summary(res)

res.id <- smacofIndDiff(perception, constraint = "identity")
res.diag <- smacofIndDiff(perception, constraint = "diagonal")
res.idio <- smacofIndDiff(perception, constraint = "idioscal")
}

\keyword{models}
