\name{EW_ger}
\alias{EW_ger}
\alias{EW_eng}
\docType{data}
\title{Work values}
\description{Intercorrelations of 13 working values for former West (first list element) and East Germany.
}

\usage{data(EW_eng)}

\format{Object of class \code{dist}
}
\details{Note that in \code{EW_ger} the labels are given in German. For smacof, the data must be converted into a dissimilarity matrix by applying the \code{sim2diss()} function to each list element.
}

\references{
ALLBUS 1991, German General Social Survey. 

Borg, I., Groenen, P. J. F., & Mair, P. (2010). Multidimensionale Skalierung. Muenchen: Hampp Verlag.

Borg, I., Groenen, P. J. F., & Mair, P. (2012). Multidimensional Scaling. New York: Springer, forthcoming.


}
\examples{
data(EW_eng)
data(EW_ger)
}
\keyword{datasets}
