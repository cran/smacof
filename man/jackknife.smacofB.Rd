\name{jackmds}
\alias{jackmds}
\alias{jackmds.smacofB}
\alias{plot.smacofJK}
\alias{print.smacofJK}


\title{SMACOF Jackknife}
\description{These methods perform a SMACOF Jackknife and plot the corresponding solution.
}
\usage{
\method{jackmds}{smacofB}(object, eps = 1e-6, itmax = 100, verbose = FALSE)

\method{plot}{smacofJK}(x, plot.dim = c(1,2), hclpar = list(c = 50, l = 70), 
col.p, col.l, plot.lines = TRUE, legend = FALSE, inset = c(-0.2, 0), cex.legend = 0.7,
main, xlab, ylab, xlim, ylim, asp = 1, ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{Object of class \code{"smacofB"}, i.e., an MDS solution from \code{smacofSym()}}
  \item{itmax}{Maximum number of iterations}
  \item{eps}{Convergence criterion}
  \item{verbose}{If \code{TRUE}, intermediate stress is printed out}
 
  \item{x}{Object of class \code{"smacofJK"}}
  \item{plot.dim}{Vector with dimensions to be plotted.}
  \item{hclpar}{Chroma and luminance to be used for HCL colors (further details see \code{\link{rainbow_hcl}})}
  \item{col.p}{Point color. If omitted, hcl colors will be used; if specified, the corresponding (single) color will be used for plotting.}
  \item{col.l}{Line color. If omitted, hcl colors will be used; if specified, the corresponding (single) color will be used for plotting.}
  \item{plot.lines}{If \code{TRUE}, the Jackknife configurations are plotted are connected with their centroid.}
  \item{legend}{If \code{TRUE}, the centroid labels are added as legend.}
  \item{inset}{Inset distance from the margins as a fraction of the plot region when legend is placed by keyword.}
  \item{cex.legend}{Character expansion factor for legend.}
  \item{main}{Plot title.}
  \item{xlab}{Label of x-axis.}
  \item{ylab}{Label of y-axis.}
  \item{xlim}{Scale x-axis.}
  \item{ylim}{Scale y-axis.}
  \item{asp}{Aspect ratio.}
  \item{\dots}{Further plot arguments passed: see \code{\link[graphics:plot.default]{plot}} for detailed information.}
}

\details{In order to examine the stability solution of an MDS, a Jackknife on the configurations can be performed (see de Leeuw & Meulman, 1986) and plotted. The plot shows the jackknife configurations which are connected to their centroid. In addition, the original smacof configuration (transformed through Procrustes) is plotted. The Jackknife function itself returns also a stability measure (as ratio of between and total variance), a measure for cross validity, and the dispersion around the original smacof solution. 
}

\value{   
  \item{smacof.conf}{SMACOF configurations}
  \item{jackknife.conf}{An array of n-1 configuration matrices for each Jackknife MDS solution}
  \item{comparison.conf}{Centroid Jackknife configurations (comparison matrix)}
  \item{stab}{Stability measure}
  \item{cross}{Cross validity}
  \item{disp}{Dispersion}
  \item{loss}{Value of the loss function}
  \item{ndim}{Number of dimensions}
  \item{call}{Model call}
  \item{niter}{Number of iterations}
  \item{nobj}{Number of objects}
}

\references{De Leeuw, J., & Meulman, J. (1986). A special jackknife for multidimensional scaling. Journal of Classification, 3, 97-112.  
}
\author{Jan de Leeuw and Patrick Mair}

\seealso{\code{\link{bootmds}}}

\examples{

## symmetric smacof
data <- na.omit(PVQ40[,1:5])
diss <- dist(t(data))   ## Euclidean distances 
fit <- mds(diss) 
res.jk <- jackmds(fit)

plot(res.jk, col.p = "black", col.l = "gray")
plot(res.jk, hclpar = list(c = 80, l = 40))
plot(res.jk, hclpar = list(c = 80, l = 40), plot.lines = FALSE)

}

\keyword{ models }
\keyword{ hplot }
