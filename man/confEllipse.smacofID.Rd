\name{confEllipse}
\alias{confEllipse}
\alias{confEllipse.smacofB}
\alias{confEllipse.smacofID}
\alias{plot.confell}


\title{Pseudo Confidence Ellipses}
\description{Computes pseudo-confidence ellipses for symmetric and individual difference MDS fits.
}
\usage{
\method{confEllipse}{smacofID}(object)

\method{plot}{confell}(x, eps = 0.05, plot.dim = c(1,2), col = 1, 
label.conf = list(label = TRUE, pos = 3, cex = 0.8), 
ell = list(lty = 1, lwd = 1, col = 1), main, xlab, ylab, xlim, ylim, 
asp = 1, type = "p", pch = 20, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{Object of class \code{"smacofB"} or \code{"smacofID"}.}
 
  \item{x}{Object of class \code{"confell"}}
  \item{eps}{Perturbation region (e.g. 0.05 means that we look at a perturbation region where stress is at most 5\% larger than the minimum we have found).}
  \item{plot.dim}{Vector with dimensions to be plotted.}
  \item{col}{Color for points.}
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
  \item{...}{Additional plotting arguments.}
}

\details{The \code{confEllipse} function normalizes the dissimilarities and performs a few more iterations to optimize the configuration and the individual diffierence weights. This result is then passed to a function that computes the stress derivatives which are the basis of the ellipses in the plot function. This function works for ratio scaled versions only. 
}

\value{   
Returns an object belonging to classes \code{"confell"}.

  \item{X}{Configuration (group stimulus space for individual difference models)}
  \item{h}{Stress derivatives}
  \item{s}{Optimized stress (raw value)}
}

\references{De Leeuw, J. (2019). Individual Differences Multidimensional Scaling. \url{http://deleeuwpdx.net/pubfolders/stability}
}



\seealso{\code{\link{plot.smacofboot}}}

\examples{
## Simple ratio MDS fit 
delta <- sim2diss(cor(PVQ40agg))
res <- mds(delta, ndim = 3)
cres <- confEllipse(res)
plot(cres, plot.dim = c(1,2))
plot(cres, plot.dim = c(1,3))
plot(cres, plot.dim = c(2,3))


## INDSCAL on Helm data
fit1 <- indscal(helm)
cfit1 <- confEllipse(fit1)
plot(cfit1, ell = list(col = "gray", lty = 2), ylim = c(-0.04, 0.04))

## IDIOSCAL on Helm data
fit2 <- idioscal(helm)
cfit2 <- confEllipse(fit2)
plot(cfit1, ell = list(col = "gray", lty = 2), ylim = c(-0.04, 0.04))

}

\keyword{hplot}
