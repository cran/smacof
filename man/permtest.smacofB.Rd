\name{permtest}
\alias{permtest}
\alias{permtest.smacofB}
\alias{plot.smacofPerm}
\alias{print.smacofPerm}


\title{SMACOF Jackknife}
\description{These methods perform a permutation test for a nonmetric symmetric SMACOF solution.
}
\usage{
\method{permtest}{smacofB}(object, nrep = 100, verbose = FALSE)

\method{plot}{smacofPerm}(x, alpha = 0.05, main, xlab, ylab, ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{Object of class \code{"smacofB"}, i.e., an MDS solution from \code{smacofSym()}}
  \item{nrep}{Number of permutations}
  \item{verbose}{If \code{TRUE}, intermediate stress is printed out}
 
  \item{x}{Object of class \code{"smacofPerm"}}
  \item{alpha}{Alpha level}
  \item{main}{Plot title.}
  \item{xlab}{Label of x-axis.}
  \item{ylab}{Label of y-axis.}
  \item{\dots}{Further plot arguments passed: see \code{\link[graphics]{plot}} for detailed information.}
}

\details{This routine permutes m dissimilarity values, where m is the number of lower diagonal elements in the corresponding dissimilarity matrix. For each sample a symmetric, nonmetric SMACOF of dimension \code{ndim} is computed and the stress values are stored in \code{stressvec}. Using the fitted stress value, the p-value is computed. Subsequently, the empirical cumulative distribution function can be plotted using the plot method. 
}

\value{   
  \item{stressvec}{Vector containing the stress values of the permutation samples}
  \item{stress.nm}{Fitted nonmetric stress}
  \item{pval}{Resulting P-values}
  \item{call}{Model call}
  \item{nrep}{Number of permutations}
  \item{nobj}{Number of objects}
}

\references{De Leeuw, J., & Stoop, I. (1984). Upper bounds for Kruskal's stress. Psychometrika, 49, 301-402.
}
\author{Jan de Leeuw and Patrick Mair}

\seealso{\code{\link{jackknife.smacofB}}, \code{\link{smacofSym}}, \code{\link{plot.smacof}}}

\examples{

## symmetric smacof
data(kinshipdelta)
res.smacof <- smacofSym(kinshipdelta, ndim = 1, metric = FALSE)
res.perm <- permtest(res.smacof, verbose = TRUE)
res.perm

plot(res.perm)
}

\keyword{models}
