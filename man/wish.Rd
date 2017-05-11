\name{wish}
\alias{wish}
\docType{data}
\title{Wish dataset}
\description{Similarity ratings for 12 countries. There were no instructions concerning the characteristics on which these similarity judgements were to be made, this was information to discover rather than to impose. 
}

\usage{data(wish)}

\format{Object of class \code{dist}
}
\details{For smacof, the data must be converted into a dissimilarity matrix (see examples). 
}

\references{
Borg, I., Groenen, P. J. F., & Mair, P. (2010). Multidimensionale Skalierung. Muenchen: Hampp Verlag.

Borg, I., Groenen, P. J. F., & Mair, P. (2012). Multidimensional Scaling. New York: Springer, forthcoming.

Wish, M. (1971). Individual differences in perceptions and preferences among nations. In C. W. King and D. Tigert (Eds.), Attitude research reaches new heights, pp. 312-328. Chicago: American Marketing Association.

}
\examples{
data(wish)
sim2diss(wish, method = max(wish))
}
\keyword{datasets}
