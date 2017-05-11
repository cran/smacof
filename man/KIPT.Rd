\name{KIPT}
\alias{KIPT}
\docType{data}
\title{Kennedy Institute Phonics Test}
\description{
Contains correlations of eight test items of the Kennedy Institute Phonics Test (KIPT), a test for reading skills. 
}
\usage{
data(KIPT)
}

\format{An 8 times 8 correlation matrix. Items: 

Nonsense word production: NP

Long vowel production: LVP

Short vowel production: SVP

Consonant cluster production: CCP 

Nonsense word recognition: NR 

Single letter production: SLP 

Consonant cluster recognition: CCR 

Initial letter recognition: ILR
}

\references{Guthrie, J. T. (1973). Models of reading and reading disability. Journal of Educational Psychology, 65, 9-18.
}

\examples{
KIPT
sim2diss(KIPT)
}
\keyword{datasets}
