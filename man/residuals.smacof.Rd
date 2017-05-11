\name{residuals.smacof}
\alias{residuals.smacof}
\alias{residuals.smacofR}
\alias{residuals.smacofID}

\title{Residuals}
\description{Computes the residuals by subtracting the configuration dissimilarites from the observed dissimilarities.
}
\usage{
\method{residuals}{smacof}(object, ...)
\method{residuals}{smacofR}(object, ...)
\method{residuals}{smacofID}(object, ...)
}
\arguments{
  \item{object}{Object of class \code{smacof}, \code{smacofR} (rectangular), or \code{smacofID} (individual differences)}
  \item{\dots}{Ignored}
}

\examples{
res <- mds(kinshipdelta)
residuals(res)

}
\keyword{methods}

