\name{summary.smacofB}
\alias{summary.smacofB}
\alias{print.smacofB}
\alias{summary.smacofR}
\alias{print.smacofR}
\alias{summary.smacofID}
\alias{print.smacofID}


\title{S3 methods for smacof}
\description{Print and summary methods for objects of class \code{smacofB}, \code{smacofR} (rectangular), and \code{smacofID} (individual differences).
}
\usage{
\method{summary}{smacofB}(object, ...)
\method{print}{smacofB}(x, ...)
\method{summary}{smacofR}(object, ...)
\method{print}{smacofR}(x, ...)
\method{summary}{smacofID}(object, ...)
\method{print}{smacofID}(x, ...)
}
\arguments{
  \item{object}{Object of class \code{smacofB}, \code{smacofR}, \code{smacofID}}
  \item{x}{Object of class \code{smacofB}, \code{smacofR}, \code{smacofID}}
  \item{\dots}{Ignored}
}

\examples{
data(kinshipdelta)
res <- smacofSym(kinshipdelta)
res
summary(res)


}
\keyword{methods}

