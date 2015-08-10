\name{sim2diss}
\alias{sim2diss}

\title{Converts similarites to dissimilarities}
\description{Utility function for converting similarities into dissimilarities. Different methods are provided.
}
\usage{
sim2diss(similmat, method = "corr", to.dist = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{similmat}{Similarity matrix (not necessarily symmetric, nor quadratic)}
  \item{method}{Various methods for converting similarities into dissimilarities: \code{"corr"}, \code{"neglog"}, \code{"counts"}, or an integer value (see details)} 
  \item{to.dist}{If \code{TRUE}, object of class dist is produced}
}
\details{We provide the following methods for converting similarities S into dissimilarities D: \code{"corr"} is suited for correlation matrices and takes \eqn{D = \sqrt{1-S}}. \code{"neglog"} takes the negative logarithm in terms of \eqn{-\log(S)}. Having frequencies, \code{"counts"} is appropriate which does \eqn{-\log((S[i,j]*S[j,i)/(S[i,i]*S[j,j]))}. The user can specify also an integer value v. In this case \code{sim2diss()} computes \eqn{v-S}. 
}
\value{Returns dissimiarities either as matrix or as dist object.
}

\examples{

## Converting Ekman data (similarities) into dissimilarities by subtraction from 1
data(ekman)
ekman.diss <- sim2diss(ekman, method = 1)
res <- smacofSym(ekman.diss)
}

\keyword{models}
