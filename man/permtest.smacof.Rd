\name{permtest}
\alias{permtest}
\alias{permtest.smacof}
\alias{permtest.smacofR}
\alias{plot.smacofPerm}
\alias{print.smacofPerm}


\title{SMACOF Permutation}
\description{These methods perform a permutation test for a symmetric or an unfolding SMACOF model.
}
\usage{
\method{permtest}{smacof}(object, data,  method.dat = "pearson", nrep = 100, verbose = TRUE, ...)
\method{permtest}{smacofR}(object, data = NULL,  method.dat = "rows", nrep = 100, verbose = TRUE, ...)
\method{plot}{smacofPerm}(x, alpha = 0.05, main, xlab, ylab, ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{Object of class \code{"smacofB"}, i.e., an MDS solution from \code{smacofSym()}}
  \item{data}{Optional argument; if provided permutations are performed on the data matrix (see details; ignored for unfolding models)}
  \item{method.dat}{If data are provided, this must be one of \code{"pearson"}, \code{"spearman"}, \code{"kendall"}, \code{"euclidean"}, \code{"maximum"}, \code{"manhattan"}, \code{"canberra"}, \code{"binary"}. For unfolding models it is either \code{"full"} for full permutations or \code{"rows"} for permutations within rows.} 
  \item{nrep}{Number of permutations}
  \item{verbose}{If \code{TRUE}, permutation index is printed out}
     
  \item{x}{Object of class \code{"smacofPerm"}}
  \item{alpha}{Alpha level}
  \item{main}{Plot title.}
  \item{xlab}{Label of x-axis.}
  \item{ylab}{Label of y-axis.}
  \item{\dots}{additional plot arguments for plot function; additional arguments to be passed to \code{sim2diss} in permutation functions.}
}

\details{This routine permutes m dissimilarity values, where m is the number of lower diagonal elements in the corresponding dissimilarity matrix. For each sample a symmetric, nonmetric SMACOF of dimension \code{ndim} is computed and the stress values are stored in \code{stressvec}. Using the fitted stress value, the p-value is computed. Subsequently, the empirical cumulative distribution function can be plotted using the plot method. 

If the MDS fit provided on derived proximities of a data matrix, this matrix can be passed to the \code{permtest} function. Consequently, the data 
matrix is subject to permutations. The proximity measure used for MDS fit has to match the one used for the permutation test. If a correlation 
similarity is provided, it is converted internally into a dissimilarity using \code{sim2diss} with corresponding arguments passed to the
\code{...} argument. 

}

\value{   
  \item{stressvec}{Vector containing the stress values of the permutation samples}
  \item{stress.obs}{Stress (observed sample)}
  \item{pval}{Resulting p-value}
  \item{call}{Model call}
  \item{nrep}{Number of permutations}
  \item{nobj}{Number of objects}
}


\author{Patrick Mair and Ingwer Borg}

\seealso{\code{\link{jackmds}}, \code{\link{bootmds}}}

\examples{

## permuting the dissimilarity matrix (full)
data(kinshipdelta)
fitkin <- mds(kinshipdelta, ndim = 2, type = "interval")
set.seed(222)
res.perm <- permtest(fitkin)
res.perm
plot(res.perm)

## permuting the data matrix
GOPdtm[GOPdtm > 1] <- 1     ## use binary version
diss1 <- dist(t(GOPdtm[,1:10]), method = "binary")  ## Jaccard distance
fitgop1 <- mds(diss1, type = "ordinal")
fitgop1
set.seed(123)
permtest(fitgop1, GOPdtm[,1:10], nrep = 10, method.dat = "binary")

rmat <- cor(GOPdtm[,1:10], method = "kendall")  ## Kendall correlation
diss2 <- sim2diss(rmat, method = 1)
fitgop2 <- mds(diss2, type = "ordinal")
fitgop2
set.seed(123)
permtest(fitgop2, GOPdtm[,1:10], nrep = 10, method.dat = "kendall", method = 1)

## unfolding permutation
data(breakfast)
res.unfolding <- unfolding(breakfast, ndim = 2)
set.seed(123)
permtest(res.unfolding, nrep = 20, method.dat = "rows")
}

\keyword{models}
