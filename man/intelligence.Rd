\name{intelligence}
\alias{intelligence}
\docType{data}
\title{Intelligence Tests}
\description{Contains intercorrelations of 8 intelligence tests, together with two facets. 
In addition, a hypothesized restriction matrix for the intercorrelations is provided. 
The proximities for items with the same structuples, such as p(NA1,NA2) and p(GA1,GA3), all are set to the value 
5. The proximities that correspond to the immediate neighborhood relations
are set to the value 4, since none of these distances should be larger 
than any distance between definitionally equivalent items. Finally, the large 
distances between the groups NI, GA and the groups NA, GI are set to 3.
The intelligence tests are coded on the follwing facets: format (N = numerical, G = geometrical) and 
requirement (A = application, I - inference). 
}
\usage{data(intelligence)}
\format{Data frame of 8 intelligence tests: facets, intercorrelations, and restrictions

Test: Test number

Languange: numerical, geometrical

Requirement: application, inference

T1-T8: intercorrelations

R1-R8: restrictions

}

\references{Borg, I., & Groenen, P. J. F. (2005). Modern Multidimensional Scaling (2nd ed.). Springer.
}
\examples{
idiss <- sim2diss(intelligence[,paste0("T", 1:8)])
fit <- mds(idiss)
plot(fit)

}
\keyword{datasets}
