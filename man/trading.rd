\name{trading}
\alias{trading}
\docType{data}
\title{Trading data}
\description{Data from the New Geographical Digest (1986), analysed in Cox and Cox (2001),
on which countries traded with other countries. For 20 countries the main trading partners 
are dichotomously scored (1 = trade performed, 0 = trade not performed). 
Based on this dichotomous matrix the dissimilarities are computed using the Jaccard coefficient.
}
\usage{data(trading)}
\format{Object of class \code{"dist"} with dissimilarities of the following countries:

Arge: Argentina

Aust: Australia

Braz: Brazil

Cana: Canada

Chin: China

Czec: Czechoslovakia 

Egyp: Egypt

E.Ge: East Germany

Fran: France

Hung: Hungary

Indi: India

Ital: Italy

Japa: Japan

N.Ze: New Zealand

Pola: Poland

Swed: Sweden

USA

USSR: Soviet Union

U.K.: United Kingdom

W.Ge: West Germany

}

\references{Cox, T.F., Cox, M.A.A. (1991). Multidimensional scaling on a sphere. Communications in Statistics: Theory and Methods, 20, 2943-2953.
}
\examples{
data(trading)

}
\keyword{datasets}
