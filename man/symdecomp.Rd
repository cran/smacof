\name{symdecomp}
\alias{symdecomp}

\title{Proximity Matrix Decomposition}
\description{Additive decomposition of an asymmetric, square proximity matrix into a symmetric matrix and an skew-symmetric matrix
}

\usage{
symdecomp(P)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{P}{Square matrix}
}

\details{
Performs the decompositon \eqn{P = M + N} (M and N are orthogonal).
}
\value{
Returns the following matrices: 
  \item{M}{symmetric component}
  \item{N}{skew-symmetric component}

}

\references{Borg, I., & Groenen, P. J. F. (2005). Modern Multidimensional Scaling (2nd ed.). Springer.
}

\examples{
P <- matrix(c(92,5,4,8,4,84,38,62,6,37,87,17,13,31,17,88), ncol = 4)
symdecomp(P)
}

\keyword{utilities}
