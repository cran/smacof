\name{stardist}
\alias{stardist}
\docType{data}
\title{Distances among stars in zodiac signs}
\description{A distance matrix for the 10 brightest stars in each of the 12 zodiac signs was computed. Astronomers measure the projected positions of objects on the celestial sphere in two angles, i.e. right ascension \eqn{\alpha} and declination \eqn{\delta}. For every zodiac sign, the projected distances on the sky between individual stars \eqn{S_{i}} and \eqn{S_{j}} have been calculated in decimal degrees by means of the Pythagorean theorem
\deqn{d_{i,j} = \sqrt {\left( {\alpha _i - \alpha _j } \right)2 + \left( {\delta _i - \delta _j } \right)2 }} assuming planar geometry. Since the zodiac signs are relatively small compared to the whole celestial sphere and the computation is only done for illustrative purposes, such a simplified assumption is appropriate.
}
\usage{data(stardist)}
\format{A dist object contaning the star distances.
}

\note{Thanks to Paul Eigenthaler, Department of Astronomy, University of Vienna for calculating the distances.}
\examples{
data(stardist)

}
\keyword{datasets}
