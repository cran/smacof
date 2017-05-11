\name{stress0}
\alias{stress0}


\title{Zero-Iterations Stress}
\description{
Computes the stress for 0 iterations based on a starting configuration provided by the user.
}
\usage{
stress0(delta, init, type = c("interval", "ratio", "ordinal", "mspline"), 
        weightmat = NULL, ties = "primary", spline.degree = 2, spline.intKnots = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{delta}{Either a symmetric dissimilarity matrix or an object of class \code{"dist"}}
  \item{init}{An initial configuration provided by the user}
  \item{weightmat}{Optional matrix with dissimilarity weights}
  \item{type}{MDS type: \code{"interval"}, \code{"ratio"}, \code{"ordinal"} (nonmetric MDS), or \code{"mspline"}}
  \item{ties}{Tie specification (ordinal MDS only): \code{"primary"}, \code{"secondary"}, or \code{"tertiary"}}
  \item{spline.degree}{Degree of the spline for \code{"mspline"} MDS type}
  \item{spline.intKnots}{Number of interior knots of the spline for \code{"mspline"} MDS type}
}
\details{Computes stress-1 for a particular starting configuration the user needs to provide. It can also be helpful if the user wants to move some points in a particular configuration such that it fits some theoretical expectations.
}

\value{Stress-1 value
}

\seealso{\code{\link{mds}}
}
\examples{
## rectangle starting solution
rect_constr
stress0(rectangles, init = rect_constr)

## torgerson starting solution
tstart <- torgerson(rectangles)
stress0(rectangles, init = tstart)
}

\keyword{utilities}
