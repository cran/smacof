\name{bootmds}
\alias{bootmds}
\alias{bootmds.smacofB}
\alias{plot.smacofboot}


\title{SMACOF Bootstrap}
\description{Performs a bootstrap on a SMACOF solution. It works for derived dissimilarities only. The original data matrix needs to 
be provided, as well as the type of dissimilarity measure used to compute the input dissimilarities.
}
\usage{
\method{bootmds}{smacofB}(object, data,  method.dat = "pearson", nrep = 100, 
alpha = 0.05, verbose = FALSE, ...)

\method{plot}{smacofboot}(x, plot.dim = c(1,2), col = 1, 
label.conf = list(label = TRUE, pos = 3, cex = 0.8), 
ell = list(lty = 1, lwd = 1, col = 1), main, xlab, ylab, xlim, ylim, asp = 1, 
type = "p", pch = 20, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{object}{Object of class \code{"smacofB"}, i.e., an MDS solution from \code{mds()}.}
  \item{data}{Initial data (before dissimilarity computation).}
  \item{method.dat}{Dissimilarity computation used as MDS input. This must be one of \code{"pearson"}, \code{"spearman"}, \code{"kendall"}, \code{"euclidean"}, \code{"maximum"}, \code{"manhattan"}, \code{"canberra"}, \code{"binary"}. For unfolding models it is either \code{"full"} for full permutations or \code{"rows"} for permutations within rows.} 
  \item{nrep}{Number of bootstrap replications.}
  \item{alpha}{Alpha level for confidence ellipsoids.}
  \item{verbose}{If \code{TRUE}, bootstrap index is printed out.}
  \item{...}{Additional arguments needed for dissimilarity computation as specified in \code{sim2diss()}.}
 
  \item{x}{Object of class \code{"smacofboot"}}
  \item{plot.dim}{Vector with dimensions to be plotted.}
  \item{col}{Color for points/ellipses.}
  \item{label.conf}{List with arguments for plotting the labels of the configurations in a configuration plot (logical value whether to plot labels or not, label position). If \code{pos = 5} labels are placed away from the nearest point.}
  \item{ell}{List with arguments for plotting ellipses: line type, line width, color.}
  \item{main}{Plot title.}
  \item{xlab}{Label of x-axis.}
  \item{ylab}{Label of y-axis.}
  \item{xlim}{Scale x-axis.}
  \item{ylim}{Scale y-axis.}
  \item{asp}{Aspect ratio.}
  \item{pch}{Plotting symbol for object point.}
  \item{type}{Type of plot.}
}

\details{In order to examine the stability solution of an MDS, a bootstrap on the raw data can be performed. This results in confidence ellipses in the configuration plot. The ellipses are returned as list which allows users to produce (and further customize) the plot by hand. 
}

\value{   
  \item{cov}{Covariances for ellipse computation}
  \item{bootconf}{Configurations bootstrap samples}
  \item{stressvec}{Bootstrap stress values}
  \item{bootci}{Stress bootstrap percentile confidence interval}
}

\references{Jacoby, W. G., & Armstrong, D. A. (2014). Bootstrap confidence regions for multidimensional scaling solutions. American Journal of Political Science, 58, 264-278. 
}


\seealso{\code{\link{jackknife}}}

\examples{

## Example using Euclidean distances
data <- na.omit(PVQ40[,1:5])
diss <- dist(t(data))   ## Euclidean distances 
fit <- mds(diss)        ## 2D interval MDS

set.seed(123)
resboot <- bootmds(fit, data, method.dat = "euclidean", nrep = 50)
resboot
plot(resboot)

## Example using Pearson correlations 
sim <- cor(data)
diss <- sim2diss(sim, method = 1)  ## subtract from 1 (method needs to be passed to bootmds)
fit <- mds(diss, type = "ratio", ndim = 3)        ## 3D ratio MDS

set.seed(123)
resboot <- bootmds(fit, data, method.dat = "pearson", nrep = 50, alpha = 0.1, method = 1)
resboot
## plot 1st against 3rd dimension
ellipses <- plot(resboot, plot.dim = c(1,3), ell = list(lty = 2, col = "gray", lwd = 0.8))
str(ellipses)  ## list of ellipse coordinates for each object
}

\keyword{ models }
\keyword{ hplot }
