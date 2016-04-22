\name{dissWeights}
\alias{dissWeights}

\title{Create Weights for Uniform Weighted Distribution}
\description{
Compute weights as a function of the dissimilarities.
}
\usage{
dissWeights(delta, type = c("unif", "knn", "power", "unifpower"), 
           k = NULL, power = 0)
}
\arguments{
  \item{delta}{Either a symmetric dissimilarity matrix or an object of class \code{"dist"}}
  \item{type}{One of \code{"unif"} (default), \code{"knn"}, \code{"power"}, \code{"unifpower"}. See details for a description of the various options.}
  \item{k}{The number of smallest dissimilarities per row for which the weights need to be set to 1. The default \code{k = NULL} makes \code{k} to be set to use the 25\% smallest dissimilarities per row.}
  \item{power}{power to which the dissimilarities need to be raised as weights. Default is 0, so that all weights are 1.}
}
\details{
The weights are computed as a function of the dissimilarities depending on \code{type}.

- \code{"unif"} Compute weights such that the weighted empirical distribution (histogram) of the dissimilarities is uniform. Particularly if the number of objects is large, the dissimilarities that occuur most often will start to dominate the solution. This option de-emphasizes often occuring dissimilarities such that the weighted empirical distribution (the weighted histogram) becomes approximately uniform.

- \code{"knn"} Per row of the dissimilarity matrix the \code{k} smallest dissimilarities obtain a weight of 1 and the others a 0.

- \code{"power"} The weights are set to the \code{delta^power}. If \code{power} is small (e.g., \code{power = -5}) then the smaller dissimilarities will be better fitted. If \code{power} is large (e.g., \code{power = 5}) then the larger dissimilarities will be better fitted.  

- \code{"unifpower"} First weights are determined by the \code{"unif"} option and then multiplied by the weights obtained by the \code{"power"} option. If the dissimilarity matrix is large, then this option is needed to see an effect of the \code{"power"} option on the MDS solution.  

}
\value{
  \item{weightmat}{the weight matrix}
}
\author{Patrick Groenen}

\examples{

## mds solution for kinship data with uniform weights
res <- mds(kinshipdelta, weightmat = dissWeights(kinshipdelta, type = "unif"))
par(mfrow = c(2,2))
plot(res, main = "uniform weights")
plot(res, plot.type = "Shepard")
plot(res, plot.type = "histogram")

## mds solution for kinship data with knn weights
res <- mds(kinshipdelta, weightmat = dissWeights(kinshipdelta, type = "knn", k = 5))
par(mfrow = c(1,2))
plot(res, main = "knn weights with k=5")
plot(res, plot.type = "Shepard")

## mds solution for kinship data with power weights emphasizing large dissimilarities
res <- mds(kinshipdelta, weightmat = dissWeights(kinshipdelta, type = "power", power = 5))
par(mfrow = c(2,2))
plot(res, main = "Power = 5 weights")
plot(res, plot.type = "Shepard")
plot(res, plot.type = "histogram")

## mds solution for kinship data with power weights emphasizing small dissimilarities
res <- mds(kinshipdelta, weightmat = dissWeights(kinshipdelta, type = "power", power = -5))
par(mfrow = c(2,2))
plot(res, main = "Power = -5 weights")
plot(res, plot.type = "Shepard")
plot(res, plot.type = "histogram")

## mds solution for kinship data with power weights emphasizing large dissimilarities 
## while correcting for nonuniform dissimilarities
res <- mds(kinshipdelta, weightmat = dissWeights(kinshipdelta, type = "unifpower", power = 5))
par(mfrow = c(2,2))
plot(res, main = "Uniform power = 5 weights")
plot(res, plot.type = "Shepard")
plot(res, plot.type = "histogram")

## mds solution for kinship data with power weights emphasizing small dissimilarities
## while correcting for nonuniform dissimilarities
res <- mds(kinshipdelta, weightmat = dissWeights(kinshipdelta, type = "unifpower", power = -5))
par(mfrow = c(2,2))
plot(res, main = "Uniform power = -5 weights")
plot(res, plot.type = "Shepard")
plot(res, plot.type = "histogram")
}
\keyword{weights}
