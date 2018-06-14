\name{inverseMDS}
\alias{inverseMDS}

\title{Inverse MDS}
\description{
Inverse MDS finds the set of all dissimilarity matrices, whose MDS solution is given as input configuration. 
}
\usage{
inverseMDS(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Configuration matrix}
}

\value{
  \item{confdist}{Euclidean distance matrix}
  \item{dissmat}{List of dissimilarity matrices}
}


\details{
One solution to inverse MDS is simply the Euclidean distances between the points in the low-dimensional space. There will typically be many additional solutions, i.e., many high-dimensional dissimilarity arrangements optimally represented by a given low-dimensional configuration. 
}

\references{
De Leeuw, J., & Groenen P. J. F. (1997). Inverse multidimensional scaling. Journal of Classification, 12, 3-21.

De Leeuw, J., Groenen P. J. F., & Mair, P. (2016). On Inverse Multidimensional Scaling. \url{http://rpubs.com/deleeuw/148292}.
}

\seealso{\code{\link{mds}}}

\examples{
## Inverse MDS on theoretical rectangle configuration
ifit <- inverseMDS(rect_constr)   
str(ifit)
}

\keyword{models}
