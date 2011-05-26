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
\details{If the constraint is \code{NULL}, INDSCAL is performed with identity restriction on the weight matrix. 
IDIOSCAL can be computed using the \code{"idioscal"} argument. The weight matrices are then unconstrained. Addtional weight restrictions can be imposed with \code{"identity"} which restricts the configurations across individuals/replications/ways to be equal, or \code{"diagonal"} which restricts only the configuration weight matrices to be diagonal. 
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
  \item{sps}{Stress per subject (matrix)}
  \item{ndim}{Number of dimensions}
  \item{model}{Type of smacof model}
  \item{niter}{Number of iterations}
  \item{nobj}{Number of objects}

}
\references{de Leeuw, J., & Mair, P. (2009). Multidimensional scaling using majorization: 
The R package smacof. Journal of Statistical Software, 31(3), 1-30, \url{http://www.jstatsoft.org/v31/i03/} 
}
\author{Jan de Leeuw and Patrick Mair}

\seealso{\code{\link{smacofConstraint}}, \code{\link{smacofSym}}, \code{\link{smacofRect}}, \code{\link{smacofSphere.primal}}, \code{\link{smacofSphere.dual}}}

\examples{
data(perception)
res <- smacofIndDiff(perception)  ## INDSCAL
res
summary(res)

res.id <- smacofIndDiff(perception, constraint = "identity")     ## identity restricted weights
res.id$cweights
plot(res.id)
plot(res.id, type = "p", pch = 25, col = 4, label.conf = list(label = TRUE, pos = 3, col = 4))


res.diag <- smacofIndDiff(perception, constraint = "diagonal")  ## diagonally restricted weights
res.diag$cweights
plot(res.diag)
plot(res.diag, type = "p", pch = 25, col = 4, label.conf = list(label = TRUE, pos = 3, col = 4))

res.idio <- smacofIndDiff(perception, constraint = "idioscal")  ## IDIOSCAL
Wk <- res.idio$cweights
G <- res.idio$gspace
G%*%Wk[[1]]  ## configurations subject 1 (cf. res.idio$conf[[1]])
G%*%Wk[[2]]  ## configurations subject 2 (cf. res.idio$conf[[3]])


data(helm)
res.helm <- smacofIndDiff(helm, constraint = "diagonal")
plot(res.helm, plot.type = "confplot")
barplot(sort(res.helm$sps, decreasing = TRUE), main = "Stress per Subject", cex.names = 0.8) 
plot(res.helm, plot.type = "bubbleplot")
plot(res.helm, plot.type = "stressplot")
plot(res.helm, plot.type = "Shepard")
plot(res.helm, plot.type = "resplot")

}

\keyword{models}
