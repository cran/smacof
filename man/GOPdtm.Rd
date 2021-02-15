\name{GOPdtm}
\alias{GOPdtm}
\docType{data}
\title{Republican Statements}
\description{Document-term matrix based on statements by Republican voters. 
}
\usage{data(GOPdtm)}
\format{Document-term matrix with statements in the rows and terms (keywords) in the columns
}
\details{This dataset emerges from statements of Republican voters scraped from the official 
GOP website. They were asked to complete the sentence "I am a Republican because ...". We have 
selected the 37 most frequent words and created a document-term matrix.
}

\references{air, P., Rusch, T. & Hornik, K. (2014). The Grand Old Party - A party of values? SpringerPlus, 3(697), \url{https://springerplus.springeropen.com/articles/10.1186/2193-1801-3-697} 
}

\examples{
data(GOPdtm)
GOPdtm
}
\keyword{datasets}
