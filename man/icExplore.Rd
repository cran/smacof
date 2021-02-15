\name{icExplore}
\alias{icExplore}
\alias{plot.icexplore}

\title{Exploring Initial Configurations}
\description{
Allows to user to explore the effect of various random starting configurations when fitting an MDS model.
}
\usage{
icExplore(delta, nrep = 100, returnfit = FALSE, ndim = 2, 
type = c("ratio", "interval", "ordinal","mspline"), weightmat = NULL, ties = "primary",	
verbose = FALSE, relax = FALSE, modulus = 1, itmax = 1000, eps = 1e-6, 
spline.degree = 2, spline.intKnots = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{delta}{Either a symmetric dissimilarity matrix or an object of class \code{"dist"}}
  \item{nrep}{Number of initial random configurations}
  \item{returnfit}{If \code{TRUE} all fitted models are returned.}
  \item{ndim}{Number of dimensions}
  \item{weightmat}{Optional matrix with dissimilarity weights}
  \item{type}{MDS type: \code{"interval"}, \code{"ratio"}, \code{"ordinal"} (nonmetric MDS), or \code{"mspline"}}
  \item{ties}{Tie specification (ordinal MDS only): \code{"primary"}, \code{"secondary"}, or \code{"tertiary"}}
  \item{verbose}{If \code{TRUE}, replication number is printed}
  \item{relax}{If \code{TRUE}, block relaxation is used for majorization}
  \item{modulus}{Number of smacof iterations per monotone regression call}
  \item{itmax}{Maximum number of iterations}
  \item{eps}{Convergence criterion}
  \item{spline.degree}{Degree of the spline for \code{"mspline"} MDS type}
  \item{spline.intKnots}{Number of interior knots of the spline for \code{"mspline"} MDS type}
}

\details{
This function generates a large set of MDS solutions using random initial configurations, matches them all by Procrustean fittings, computes the inter-correlations of their point coordinates, and finally runs an interval MDS of these inter-correlations. It can be used to explore local minima.

In the plot function the number reflects the index of corresponding MDS fit, the size reflects the stress value: the larger the font, the larger the stress (i.e., the worse the solution). The size is associated with a corresponding color shading (the smaller the size the darker the color). 
}

\value{
  \item{mdsfit}{Fitted MDS objects (\code{NULL} if \code{returnfit = FALSE})}
  \item{conf}{Configuration based on multiple random starts}
  \item{stressvec}{Vector with stress values}
}

\references{
Borg, I. and Mair, P. (2017). The choice of initial configurations in multidimensional scaling: local
minima, fit, and interpretability. Austrian Journal of Statistics, 46, 19-32. \doi{10.17713/ajs.v46i2.561}
}


\seealso{\code{\link{mds}}}
\examples{

## simple example with 20 random starts
diss <- sim2diss(wish, method = 7)
set.seed(123)
res <- icExplore(diss, type = "ordinal", nrep = 20, returnfit = TRUE)
res
plot(res) 

res$mdsfit[[14]]  ## bad fitting solution
res$mdsfit[[3]]   ## better fitting solution
}

\keyword{utilities}
