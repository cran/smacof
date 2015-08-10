\name{uniscale}
\alias{uniscale}
\alias{print.uniscale}
\alias{plot.uniscale}

\title{Unidimensional Scaling}
\description{
Simple implementation where all dissimilarity permutations are subject to a 1D MDS fit and the one which leads to a minimal stress values is returned.
}
\usage{
uniscale(delta, weightmat = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{delta}{Either a symmetric dissimilarity matrix or an object of class \code{"dist"}}
  \item{weightmat}{Optional matrix with dissimilarity weights}
}

\value{
  \item{delta}{Observed dissimilarities, not normalized}
  \item{confdiss}{Configuration distances}
  \item{conf}{Vector with fitted configurations}
  \item{stress}{Stress-1 value}
  \item{weightmat}{Weight matrix}
  \item{nobj}{Number of objects}
  \item{npermtot}{Total number of permutations (factorial)}
  \item{npermscale}{Number of accepted permutations (monotonicity check)}
  
}
\references{
Mair P., De Leeuw J. (2015). Unidimensional scaling. In Wiley StatsRef: Statistics Reference Online, Wiley, New York.
}

\seealso{\code{\link{mds}}}

\examples{
## unidimensional scaling of Plato's 7 works
PlatoD <- dist(t(Plato7))
fit.uni <- uniscale(PlatoD)
fit.uni
plot(fit.uni)
}

\keyword{models}
