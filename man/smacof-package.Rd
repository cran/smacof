\name{smacof-package}
\alias{smacof-package}
\alias{smacof}
\docType{package}
\title{SMACOF for Multidimensional Scaling.
}
\description{This package provides the following approaches of multidimensional scaling (MDS) based on stress minimization by means of majorization (smacof): Simple smacof on symmetric dissimilarity matrices, smacof for rectangular matrices (unfolding models), smacof with constraints on the configuration, three-way smacof for individual differences (including constraints for idioscal, indscal, and identity), and spherical smacof (primal and dual algorithm). Each of these approaches is implemented in a metric and nonmetric manner including primary, secondary, and tertiary approaches for tie handling. 
}
\details{
\tabular{ll}{
Package: \tab smacof\cr
Type: \tab Package\cr
Version: \tab 1.4-0\cr
Date: \tab 2014-05-12\cr
License: \tab GPL\cr
}
The function for basic SMACOF on symmetric dissimilarity matrices is \code{smacofSym()}. For rectangular input matrices (unfolding model) \code{smacofRect()} is appropriate and by means of \code{smacofIndDiff()} individual difference models (three-way MDS) can be computed.

}
\author{
Jan de Leeuw, Patrick Mair

Maintainer: Jan de Leeuw <deleeuw@stat.ucla.edu>

}
\references{de Leeuw, J. \& Mair, P. (2009). Multidimensional scaling using majorization: 
The R package smacof. Journal of Statistical Software, 31(3), 1-30, \url{http://www.jstatsoft.org/v31/i03/} 
}

\keyword{ package }
\seealso{\code{\link{smacofSym}}, \code{\link{smacofRect}}, \code{\link{smacof}}, \code{\link{smacofIndDiff}}, \code{\link{smacofSphere}}
}


\examples{
data(trading)
res <- smacofSym(trading)
res
}
