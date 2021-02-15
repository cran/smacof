\name{gravity}
\alias{gravity}

\title{Gravity dissimilarities}
\description{Computes the dissimilarities using a gravity model based on co-occurrences.}

\usage{
gravity(X, lambda = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{numeric matrix}
  \item{lambda}{tuning parameter}
}

\details{The first step in this function is to compute the co-occurences. Based on the 
binarized data matrix \eqn{Y} we compute \eqn{Y'Y} which leads to the co-occurence matrix.
We then use the gravity model to compute the gravity dissimilarities. 
In order to give more (or less) structure to the MDS solution, the tuning parameter (which 
defines a power transformation) can be increased (or decreased). In addition, 
a weight matrix is created that sets cells with no co-occurences to 0. The corresponding weight matrix for blanking out these cells is established automatically in \code{mds()}. 
}

\value{
  \item{gravdiss}{Gravity dissimilarities}
  \item{weightmat}{Weight matrix for subsequent smacof computation}
  \item{co.occ}{Matrix with co-occurences}
}
\references{Mair, P., Rusch, T. & Hornik, K. (2014). The Grand Old Party - A party of values?
SpringerPlus, 3(697), \url{https://springerplus.springeropen.com/articles/10.1186/2193-1801-3-697} 
}
\author{Patrick Mair}

\seealso{\code{\link{mds}}}

\examples{
data(GOPdtm)
gravD <- gravity(GOPdtm, lambda = 2)
res <- mds(gravD$gravdiss)  
res$weightmat  ## NA's were blanked out when fitting the model
plot(res)
}

\keyword{models}
