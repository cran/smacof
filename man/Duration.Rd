\name{Duration}
\alias{Duration}
\alias{DurationRaw}
\docType{data}
\title{Facets of the subjective duration of imagined situations}
\description{The \code{DurationRaw} dataset contains the duration rating of 76 subjects on 24 situations. 
Subjects were asked to rate the duration on a 7 point scale (1 ... substantially shorter, 
7 ... substantially longer). The \code{Duration} data file contains the corresponding 
correlations between the 24 situations including some information about the facets. 
}
\usage{
data(Duration)
data(DurationRaw)
}
\format{Data frame 24 correlations based on duration ratings:

S1-S24: situation

F1: pleasant (1), neutral (2), unpleasant (3)

F2: variable (1), monotonous (2)

F3: difficult (1), easy (2)

F3: many (1), few (2)

structuple: the facet structure written as a tuple
}

\references{Borg, I., & Groenen, P. J. F. (2005). Modern Multidimensional Scaling (2nd ed.). Springer.
}
\examples{
ddiss <- sim2diss(Duration[,paste0("S", 1:24)])
fit <- mds(ddiss, type = "ordinal", ndim = 4)
plot(fit)

}
\keyword{datasets}
