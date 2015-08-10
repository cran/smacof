\name{torgerson}
\alias{torgerson}

\title{Torgerson Scaling}
\description{Classical MDS aka Torgerson Scaling
}

\usage{
torgerson(delta, p)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{delta}{Dissimilarity matrix}
  \item{p}{Number of dimensions}
}


\value{Returns an \eqn{n \times p} matrix of configurations
}

\references{Borg, I., & Groenen, P. J. F. (2005). Modern Multidimensional Scaling (2nd ed.). Springer.
}

\examples{
fit <- torgerson(Guerry)
}

\keyword{models}
