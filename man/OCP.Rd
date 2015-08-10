\name{OCP}
\alias{OCP}
\docType{data}
\title{Organizational Culture Profile}
\description{Contains similarities (correlations) of 54 OCP (see O'Reilly, Chatman, and Caldwell, 1991) items. The last column contains the facet assigned by Bilsky and Jehn (2002).
}
\usage{data(OCP)}
\format{Data frame with OCP item correlations and facet:

i1-i54: OCP item correlations

facet: factor with facets
}

\references{Bilsky, W. & Jehn, K. (2002).  Organizational Culture and Individual Values: Evidence 
for a Common Structure.  In M. Myrtek (Ed.), The Person in Biological and Social Context, 
pp. 211-228. Goettingen, Germany: Hogrefe Press. 
}
\examples{
ocpD <- sim2diss(OCP[,-55])
fit <- mds(ocpD, type = "ordinal")
plot(fit)
}
\keyword{datasets}
