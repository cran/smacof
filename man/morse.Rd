\name{morse}
\alias{morse}
\alias{morse2}
\docType{data}
\title{Morse Code Confusion Data}
\description{Confusion percentages between Morse code signals. The scores are derived 
from confusion rates on 36 Morse code signals (26 for the alphabet; 10 for the numbers 
0,...,9). Each Morse code signal is a sequence of up to five 'beeps'. The beeps can be 
short (0.05 sec) or long (0.15 sec), and, when there are two or more beeps in a signal, 
 they are separated by periods of silence (0.05 sec).

Rothkopf asked 598 subjects to judge whether two signals, presented acoustically one 
after another, were the same or not. The values are the average percentages with which 
the answer 'Same!' was given in each combination of row stimulus i and column stimulus j, 
where either i or j was the first signal presented. The values are 1 minus the 
symmetrized confusion rates and are thus dissimilarities.}

\usage{
data(morse)
data(morse2)
}

\format{Similarity matrix of 36 morse codes}

\details{
The first dataset (\code{morse}) contains a symmetric version, the second dataset (\code{morse2}) the original asymmetric version.
}

\references{
Rothkopf, E. Z. (1957). A measure of stimulus similarity and errors in some 
paired-associate learning. Journal of Experimental Psychology, 53, 94-101.
}

\examples{
data(morse)

}
\keyword{datasets}
