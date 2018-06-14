\name{FaceExp}
\alias{FaceExp}
\alias{FaceScale}
\docType{data}
\title{Facial Expression Data}
\description{Dissimilarity matrix of 13 facial expressions (Abelson & Sermat, 1962). The external scales are taken from Engen et al. (1958) reflecting the following three perceptual dimensions: pleasant-unpleasant (PU), attention-rejection (AR), and tension-sleep (TS). 
}

\usage{
data(FaceExp)
data(FaceScale)
}

\format{Symmetric dissimilarity matrix and data frame with 3 perceptual dimensions}


\references{
Abelson, R. P., & Sermat, V. (1962). Multidimensional scaling of facial expressions. Journal of Experimental Psychology, 63, 546-554. 

Engen, B., Levy, N., & Schlossberg, H. (1958). The dimensional analysis of a new series of facial expressions. Journal of Experimental Psychology, 55, 454-458. 
}

\examples{
str(FaceExp)
str(FaceScale)
}
\keyword{datasets}
