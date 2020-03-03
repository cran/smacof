\name{Guttman1991}
\alias{Guttman1991}
\alias{Guttman1965}
\docType{data}
\title{Guttman's Intelligence Facets}
\description{The first dataset (\code{Guttman1991}) contains similarities and facets for Guttman's 3D cylindrical intelligence structure as published in Guttman (1991). The second dataset (\code{Guttman1965}) contains similarities and structural intelligence facets from Guttman (1965). 
}

\usage{
data(Guttman1991)
data(Guttman1965)
}

\format{List with two elements: The first element contains the similarity matrix, the second element the facets labels.
}

\references{
Guttman, L. & Levy, S. (1991). Two structural laws for intelligence tests. Intelligence, 15, 79-103. 

Guttman, L. (1965). The structure of interrelations among intelligence tests. In C. W. Harris (Ed.), Proceedings of the 1964 Invitational Conference on Testing Problems (pp. 23-36). Princeton: ETS. 
}
\examples{
Guttman1991[[1]]        ## similarity matrix
Guttman1991[[2]]        ## facets

Guttman1965[[1]]        ## similarity matrix
Guttman1965[[2]]        ## facets
}
\keyword{datasets}
