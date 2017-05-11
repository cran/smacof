\name{indvalues}
\alias{indvalues}
\docType{data}
\title{Individual Psychological Values}
\description{Responses from a sample in Britain were collected varying in value measures of the Schwartz value theory. 
The instruments used was the Schwartz Value Survey (SVS).
}

\usage{data(indvalues)}
\format{Data frame with 327 persons in the rows and value dissimilarity scores in the columns:

PO: power

AC: achievement

HE: hedonism

ST: stimulation

SD: self-direction

UN: universalism

BE: benevolence 

TR: tradition

CO: conformity  

SE: security
}

\details{
The data were centered (row-wise) and converted from preferences into dissimilarities.
}

\references{Borg, I., Bardi, A., & Schwartz, S. H. (2017). Does the value circle exist within persons or only across persons? 
Journal of Personality, 85(2), 151-162. 
}
\examples{
str(indvalues)
}
\keyword{datasets}
