\name{smacofIndDiff}
\alias{smacofIndDiff}
\alias{indscal}
\alias{idioscal}

\title{SMACOF for Individual Differences}
\description{Performs smacof for individual differences also known as Three-Way smacof on a list of dissimilarity matrices. Various restrictions decompositions and restrictions on the weight matrix are provided. The most prominent models are 
INDSCAL and IDIOSCAL.
}

\usage{
smacofIndDiff(delta, ndim = 2, type = c("ratio", "interval", "ordinal", "mspline"), 
              constraint = c("indscal", "idioscal", "identity"),
              weightmat = NULL, init = "torgerson", ties = "primary",  
              verbose = FALSE, modulus = 1, itmax = 1000, eps = 1e-6,
              spline.degree = 2, spline.intKnots = 2)
              
indscal(delta, ndim = 2, type = c("ratio", "interval", "ordinal", "mspline"), 
        weightmat = NULL, init = "torgerson", ties = "primary", 
        verbose = FALSE, modulus = 1, itmax = 1000, eps = 1e-6,
        spline.degree = 2, spline.intKnots = 2)

idioscal(delta, ndim = 2, type = c("ratio", "interval", "ordinal", "mspline"), 
         weightmat = NULL, init = "torgerson", ties = "primary", 
         verbose = FALSE, modulus = 1, itmax = 1000, eps = 1e-6,
         spline.degree = 2, spline.intKnots = 2)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{delta}{A list of dissimilarity matrices or a list objects of class \code{dist}}
  \item{ndim}{Number of dimensions}
  \item{type}{MDS type: \code{"interval"}, \code{"ratio"}, \code{"ordinal"} (nonmetric MDS), or \code{"mspline"}}
  \item{weightmat}{Optional matrix with dissimilarity weights}
  \item{init}{Matrix with starting values for configurations (optional)}
  \item{ties}{Tie specification for non-metric MDS}
  \item{constraint}{Either \code{"indscal"}, \code{"idioscal"}, or \code{"identity"} (see details)}
  \item{verbose}{If \code{TRUE}, intermediate stress is printed out}
  \item{modulus}{Number of smacof iterations per monotone regression call}
  \item{itmax}{Maximum number of iterations}
  \item{eps}{Convergence criterion}
   \item{spline.degree}{Degree of the spline for \code{"mspline"} MDS type}
  \item{spline.intKnots}{Number of interior knots of the spline for \code{"mspline"} MDS type}
}

\details{
If the constraint is \code{"indscal"}, INDSCAL is performed with configuration weight matrices restricted to be diagonal. \code{indscal()} is a corresponding wrapper function that can be used instead of \code{smacofIndDiff()} with 
\code{"indscal"} constraints.

IDIOSCAL can be computed using the \code{"idioscal"} argument. The weight matrices are then unconstrained. 
\code{idioscal()} is a corresponding wrapper function that can be used instead of \code{smacofIndDiff()} with 
\code{"idioscal"} constraints.

Addtional weight restrictions can be imposed with \code{"identity"} which restricts the configurations across individuals/replications/ways to be equal. 
}

\value{
  \item{delta}{Observed dissimilarities}
  \item{obsdiss}{List of observed dissimilarities, normalized}
  \item{confdist}{List of configuration dissimilarities}
  \item{conf}{List of matrices of final configurations}
  \item{gspace}{Joint configuration aka group stimulus space}
  \item{cweights}{Individual weights}
  \item{stress}{Stress-1 value}
  \item{resmat}{Matrix with squared residuals}
  \item{rss}{Residual sum-of-squares}
  \item{spp}{Stress per point (in percent)}
  \item{spps}{Stress per point per subject (in percent, conditional on subject)}
  \item{sps}{Stress per subject (in percent)}
  \item{ndim}{Number of dimensions}
  \item{model}{Type of smacof model}
  \item{niter}{Number of iterations}
  \item{nobj}{Number of objects}
}

\references{
De Leeuw, J., & Mair, P. (2009). Multidimensional scaling using majorization: 
The R package smacof. Journal of Statistical Software, 31(3), 1-30, \url{http://www.jstatsoft.org/v31/i03/} 
}
\author{Jan de Leeuw and Patrick Mair}

\seealso{\code{\link{smacofConstraint}}, \code{\link{smacofSym}}, \code{\link{smacofRect}}, \code{\link{smacofSphere}}}

\examples{
## Example 1: rectangle perception data
res.diag <- indscal(perception, type = "ordinal")        ## INDSCAL
res.diag$cweights
plot(res.diag)
plot(res.diag, type = "p", pch = 25, col = 4, label.conf = list(label = TRUE, pos = 3, col = 4))

res.idio <- idioscal(perception, type = "ordinal")  ## IDIOSCAL
Wk <- res.idio$cweights
G <- res.idio$gspace
G%*%Wk[[1]]  ## configuration subject 1 (cf. res.idio$conf[[1]])
G%*%Wk[[2]]  ## configuration subject 2 (cf. res.idio$conf[[3]])

## identity restricted weights
res.id <- smacofIndDiff(perception, type = "ordinal", constraint = "identity")     
summary(res.id)
res.id$cweights
plot(res.id)
plot(res.id, type = "p", pch = 25, col = 4, label.conf = list(label = TRUE, pos = 3, col = 4))


## Example 2: Helm's color data 
res.helm <- indscal(helm, type = "interval")
plot(res.helm, plot.type = "confplot")
barplot(sort(res.helm$sps, decreasing = TRUE), main = "Stress per Subject", cex.names = 0.8) 
plot(res.helm, plot.type = "bubbleplot")
plot(res.helm, plot.type = "stressplot")
plot(res.helm, plot.type = "Shepard")

## idioscal and indscal with random starting configuration:
set.seed(123)
startconf <- matrix(rnorm(20), 10, 2)
idioscal(helm, init = startconf, type = "interval")
indscal(helm, init = startconf, type = "interval")

}

\keyword{multivariate}
