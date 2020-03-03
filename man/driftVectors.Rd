\name{driftVectors}
\alias{driftVectors}
\alias{plot.driftvec}

\title{Asymmetric MDS: Drift Vectors}
\description{
Takes an asymmetric dissimilarity matrix and decomposes it into a symmetric and a skew-symmetric part. 
Fits an MDS on the symmetric part and computes drift vectors for the skew-symmetric portion. This model makes 
it possible to see how these two components are related to each other. It is limited to two dimensions only.
}
\usage{
driftVectors(data, type = c("ratio", "interval", "ordinal","mspline"), 
             weightmat = NULL, init = "torgerson", ties = "primary",  verbose = FALSE, 
             relax = FALSE, modulus = 1, itmax = 1000, eps = 1e-6, 
             spline.degree = 2, spline.intKnots = 2)
             
\method{plot}{driftvec}(x, adjust = 1, main, xlim, ylim, 
                        xlab = "Dimension 1", ylab = "Dimension 2", pch = 20, asp = 1, 
                        col.conf = "black", col.drift = "lightgray", 
                        label.conf = list(label = TRUE, pos = 3, col = "black", 
                        cex = 0.8), ...) 

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{Asymmetric dissimilarity matrix}
    \item{weightmat}{Optional matrix with dissimilarity weights}
  \item{init}{Either \code{"torgerson"} (classical scaling starting solution), \code{"random"} (random configuration), or a user-defined matrix}
  \item{type}{MDS type: \code{"interval"}, \code{"ratio"}, \code{"ordinal"} (nonmetric MDS), or \code{"mspline"}}
  \item{ties}{Tie specification for ordinal MDS only: \code{"primary"}, \code{"secondary"}, or \code{"tertiary"}}
  \item{verbose}{If \code{TRUE}, intermediate stress is printed out}
  \item{relax}{If \code{TRUE}, block relaxation is used for majorization}
  \item{modulus}{Number of smacof iterations per monotone regression call}
  \item{itmax}{Maximum number of iterations}
  \item{eps}{Convergence criterion}
  \item{spline.degree}{Degree of the spline for \code{"mspline"} MDS type}
  \item{spline.intKnots}{Number of interior knots of the spline for \code{"mspline"} MDS type}
    
  \item{x}{Object of class \code{"driftvec"}}
  \item{adjust}{Scaling factor for drift vectors (value larger than 1 increases the length)}
  \item{main}{Plot title}
  \item{xlab}{Label of x-axis}
  \item{ylab}{Label of y-axis}
  \item{xlim}{Scale x-axis}
  \item{ylim}{Scale y-axis}
  \item{pch}{Plot symbol}
  \item{asp}{Aspect ratio}
  \item{col.conf}{Point color (MDS configurations)}
  \item{col.drift}{Color for drift vectors (arrows)}
  \item{label.conf}{Settings for plotting labels}
  \item{...}{Additional plotting arguments}
}
\details{
The skew-symmetric values are embedded into the MDS representation of the symmetrized data by drawing errors (drift vectors)
from each point \eqn{i} to each point \eqn{j} in the configuration so that these vectors correspond in length and 
direction to the values of row \eqn{i} of the skew-symmetric matrix.
}
\value{
  \item{fitsym}{MDS output for symmetric portion}
  \item{sym}{Symmetric matrix}
  \item{skewsym}{Skew-symmetric matrix}
  \item{driftcoor}{Drift vector coordinates}
  \item{stress}{Stress-1 value}
  \item{niter}{Number of iterations}
  \item{nobj}{Number of objects}
}

\references{
Borg, I., & Groenen, P. J. F. (2005). Modern Multidimensional Scaling (2nd ed.). Springer.
}

\author{Patrick Mair}

\seealso{\code{\link{smacofSym}}
}

\examples{
## simple example
P <- matrix(c(0,  4,  6, 13, 
              5,  0, 37, 21, 
              4, 38,  0, 16, 
              8, 31, 18,  0), nrow = 4, ncol = 4, byrow = TRUE) 
D <- sim2diss(P, method = 40) 
res <- driftVectors(D, type = "interval")
plot(res)
plot(res, adjust = 0.1)    ## shorten drift vectors 

## Morse code data
fit.drift <- driftVectors(morse2, type = "ordinal")
fit.drift
plot(fit.drift)
}
\keyword{multivariate}
