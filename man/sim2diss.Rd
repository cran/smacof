\name{sim2diss}
\alias{sim2diss}

\title{Converts similarites to dissimilarities}
\description{Utility function for converting similarities into dissimilarities. Different methods are provided.
}
\usage{
sim2diss(s, method = "corr", to.dist = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{s}{Similarity matrix (not necessarily symmetric, nor square)}
  \item{method}{Various methods for converting similarities into dissimilarities: \code{"corr"}, \code{"reverse"}, \code{"reciprocal"}, \code{"ranks"}, \code{"exp"}, \code{"Gaussian"}, \code{"cooccurrence"}, \code{"gravity"}, \code{"confusion"}, \code{"transition"}, \code{"membership"}, \code{"probability"}, or an integer value from which the similarity is subtracted} 
  \item{to.dist}{If \code{TRUE}, object of class dist is produced}
}
\details{The conversion formulas for the various methods can be found in the package vignette. 
}
\value{Returns dissimiarities either as matrix or as dist object.
}

\examples{

## Convert crimes data (correlations)
data(crimes)
crimeD <- sim2diss(crimes, method = "corr", to.dist = TRUE)

## Convert Wish data (similarities) by subtracting from 7
data(wish)
wishD <- sim2diss(wish, method = 7, to.dist = TRUE)

## Convert Ekman data (similarities) into dissimilarities
data(ekman)
ekmanD <- sim2diss(ekman, method = "confusion", to.dist = TRUE)

## Convert album ratings (rectangular similarities) by reversing the ratings
data(RockHard)
rockD1 <- sim2diss(RockHard[,5:18], method = "reverse")
}

\keyword{models}
